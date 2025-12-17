#include "KDTree.H"
#include "IOstreams.H"
#include <algorithm> // For std::swap

namespace Foam
{

// ------------------------------
// Constructor
// ------------------------------
KDTree::KDTree(const List<scalarField>& points)
:
    root_(),  // Initialize empty autoPtr
    points_(points),
    k_(points.empty() ? 0 : points[0].size())
{
    if (points_.empty())
    {
        return;
    }

    // Safety check for dimensionality
    forAll(points_, i)
    {
        if (points_[i].size() != k_)
        {
            FatalErrorInFunction
                << "Point " << i << " has dimension " << points_[i].size()
                << " but expected " << k_
                << abort(FatalError);
        }
    }

    // Create the master index list [0, 1, 2, ..., N-1]
    List<label> indices(points_.size());
    forAll(indices, i)
    {
        indices[i] = i;
    }

    // Build the tree using the full range
    root_ = buildTree(indices, 0, indices.size(), 0);
}


// ------------------------------
// Helper: Partition Indices (QuickSelect)
// ------------------------------
void KDTree::partitionIndices
(
    List<label>& indices,
    label start,
    label end,
    const label axis,
    const label k
)
{
    // Implementation of QuickSelect (Hoare's selection algorithm)
    // Rearranges 'indices' in [start, end) such that:
    // 1. indices[k] contains the element that would be there if sorted.
    // 2. All elements in [start, k) are <= indices[k]
    // 3. All elements in (k, end) are >= indices[k]

    label left = start;
    label right = end - 1;

    while (left < right)
    {
        const label mid = (left + right) / 2;
        
        // Pivot selection: use the middle element
        // (Swapping to right makes partition loop simpler)
        const scalar pivotVal = points_[indices[mid]][axis];
        std::swap(indices[mid], indices[right]);

        label storeIdx = left;
        for (label i = left; i < right; ++i)
        {
            if (points_[indices[i]][axis] < pivotVal)
            {
                std::swap(indices[i], indices[storeIdx]);
                storeIdx++;
            }
        }
        
        // Move pivot to its final place
        std::swap(indices[storeIdx], indices[right]);

        // Optimization: only recurse into the part that contains 'k'
        if (storeIdx == k)
        {
            return;
        }
        else if (k < storeIdx)
        {
            right = storeIdx - 1;
        }
        else
        {
            left = storeIdx + 1;
        }
    }
}


// ------------------------------
// Tree Construction (Recursive)
// ------------------------------
autoPtr<KDTree::Node> KDTree::buildTree
(
    List<label>& indices,
    const label start,
    const label end,
    const label depth
)
{
    const label count = end - start;

    if (count <= 0)
    {
        return autoPtr<Node>();
    }

    // Case 1: Leaf Node
    if (count <= leafSize_)
    {
        List<label> leafInds(count);
        for(label i = 0; i < count; ++i)
        {
            leafInds[i] = indices[start + i];
        }
        return autoPtr<Node>(new Node(leafInds));
    }

    // Case 2: Internal Node
    const label axis = depth % k_;

    // Calculate the absolute index of the median
    const label midOffset = count / 2;
    const label midIndex = start + midOffset;

    // Use QuickSelect to place the median at midIndex
    // This is O(N) instead of O(N log N) for full sort
    partitionIndices(indices, start, end, axis, midIndex);
    
    // Create internal node with the median point
    autoPtr<Node> node(new Node(indices[midIndex], axis));

    // Recursively build children
    // Left: [start, midIndex)
    node->left = buildTree(indices, start, midIndex, depth + 1);
    
    // Right: [midIndex + 1, end)
    node->right = buildTree(indices, midIndex + 1, end, depth + 1);

    return node;
}


// ------------------------------
// Distance Computation
// ------------------------------
scalar KDTree::distSqr
(
    const scalarField& a,
    const scalarField& b
) const
{
    scalar sum = 0.0;
    for (label d = 0; d < k_; ++d)
    {
        const scalar diff = a[d] - b[d];
        sum += diff*diff;
    }
    return sum;
}


// ------------------------------
// Nearest Search (Recursive)
// ------------------------------
void KDTree::nearestSearch
(
    const Node* node,
    const scalarField& query,
    label& bestIndex,
    scalar& bestDistSqr
) const
{
    if (!node) return;

    // -- Case A: Leaf Node --
    if (!node->leafIndices.empty())
    {
        forAll(node->leafIndices, i)
        {
            const label idx = node->leafIndices[i];
            const scalar d2 = distSqr(points_[idx], query);

            if (d2 < bestDistSqr)
            {
                bestDistSqr = d2;
                bestIndex = idx;
            }
        }
        return;
    }

    // -- Case B: Internal Node --
    
    // Check distance to the split point itself
    const label idx = node->pointIndex;
    const scalar d2 = distSqr(points_[idx], query);

    if (d2 < bestDistSqr)
    {
        bestDistSqr = d2;
        bestIndex = idx;
    }

    // Determine which side to search first
    const label ax = node->axis;
    const scalar diff = query[ax] - points_[idx][ax];
    
    // Identify Near and Far branches
    const autoPtr<Node>& nearBranch = (diff < 0) ? node->left : node->right;
    const autoPtr<Node>& farBranch  = (diff < 0) ? node->right : node->left;

    // 1. Search the near side
    if (nearBranch.valid())
    {
        nearestSearch(&nearBranch(), query, bestIndex, bestDistSqr);
    }

    // 2. Search the far side ONLY if the hypersphere intersects the plane
    if (farBranch.valid())
    {
        if ((diff * diff) < bestDistSqr)
        {
            nearestSearch(&farBranch(), query, bestIndex, bestDistSqr);
        }
    }
}


// ------------------------------
// Public API
// ------------------------------
label KDTree::nearest(const scalarField& query) const
{
    if (!root_.valid())
    {
        return -1;
    }

    if (query.size() != k_)
    {
        FatalErrorInFunction
            << "Query dimension " << query.size()
            << " does not match KDTree dimension " << k_
            << abort(FatalError);
    }

    label bestIndex = -1;
    scalar bestDistSqr = GREAT;

    nearestSearch(&root_(), query, bestIndex, bestDistSqr);

    return bestIndex;
}

scalar KDTree::nearestDistSqr(const scalarField& query) const
{
    if (!root_.valid())
    {
        return GREAT;
    }

    if (query.size() != k_)
    {
        FatalErrorInFunction
            << "Query dimension " << query.size()
            << " does not match KDTree dimension " << k_
            << abort(FatalError);
    }

    label bestIndex = -1;
    scalar bestDistSqr = GREAT;

    nearestSearch(&root_(), query, bestIndex, bestDistSqr);

    return bestDistSqr;
}

} // End namespace Foam