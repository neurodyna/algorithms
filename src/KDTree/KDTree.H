#ifndef KDTree_H
#define KDTree_H

#include "List.H"
#include "scalarField.H"
#include "autoPtr.H"
#include "label.H"
#include "scalar.H"
#include "pointField.H"

namespace Foam
{

class KDTree
{
private:

    // Tunable parameter: bucket size for leaf nodes
    static const label leafSize_ = 16;

    struct Node
    {
        label pointIndex;          // Point index (Internal node only)
        label axis;                // Split axis (Internal node only)
        List<label> leafIndices;   // List of points (Leaf node only)

        autoPtr<Node> left;
        autoPtr<Node> right;

        // Internal node constructor
        Node(label pI, label ax)
        :
            pointIndex(pI),
            axis(ax),
            leafIndices(),
            left(nullptr),
            right(nullptr)
        {}

        // Leaf node constructor
        Node(const List<label>& inds)
        :
            pointIndex(-1),
            axis(-1),
            leafIndices(inds),
            left(nullptr),
            right(nullptr)
        {}
    };

    // Tree Root
    autoPtr<Node> root_;

    // Reference to external points (Caller must ensure lifetime!)
    const List<scalarField>& points_;

    // Dimensionality (e.g., 3 for 3D points)
    const label k_;


    // ---- Private Helper Functions ----

    // Efficiently partitions the indices such that the element at 'k' 
    // is in its sorted position, with smaller elements to the left 
    // and larger to the right.
    void partitionIndices
    (
        List<label>& indices,
        label start,
        label end,
        const label axis,
        const label k
    );

    // Recursive builder using index ranges
    autoPtr<Node> buildTree
    (
        List<label>& indices,
        const label start,
        const label end,
        const label depth
    );

    // Squared Euclidean distance
    scalar distSqr
    (
        const scalarField& a,
        const scalarField& b
    ) const;

    // Recursive nearest neighbor search
    void nearestSearch
    (
        const Node* node,
        const scalarField& query,
        label& bestIndex,
        scalar& bestDistSqr
    ) const;

public:

    // Constructor
    explicit KDTree(const List<scalarField>& points);

    // Disable copying
    KDTree(const KDTree&) = delete;
    void operator=(const KDTree&) = delete;

    // Destructor
    ~KDTree() = default;

    // Find nearest neighbor index
    label nearest(const scalarField& query) const;

    // Find nearest neighbor squared distance
    scalar nearestDistSqr(const scalarField& query) const;
};

} // End namespace Foam

#endif