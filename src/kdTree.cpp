#include "kdTree.h"
#include <iostream>

inline void swap(std::vector<int> &v, int i, int j)
{
  int tmp = v[i];
  v[i] = v[j];
  v[j] = tmp;
}

struct sorter
{
  const Eigen::MatrixXd &points;
  int splitDim;

  sorter(const Eigen::MatrixXd &points, int splitDim) : points(points),
                                                        splitDim(splitDim) { }

  bool operator()(int i, int j)
  {
    bool result = points(i, splitDim) < points(j, splitDim);
    return result;
  }
};

int KDTreeNode::partition5(std::vector<int> &indices, int start, int end)
{
  sorter s(points, splitDim);
  std::sort(indices.begin() + start, indices.begin() + end, s);
  return (end - start) / 2 + start;
}

int KDTreeNode::getPivot(std::vector<int> &indices, int start, int end)
{
  if (end - start < 6) {
    return partition5(indices, start, end);
  }
  for (int i = start; i < end; i += 5) {
    int subEnd = i + 5;
    if (subEnd > end) subEnd = end;
    int median5 = partition5(indices, i, subEnd);
    swap(indices, median5, start + (i - start) / 5);
  }
  return findMedian(indices, start, start + (end - start + 4) / 5);
}

int KDTreeNode::findMedian(std::vector<int> &indices, int start, int end)
{
  if (end - start == 1) return start;
  int n = (end - start) / 2 + start;
  while (true) {
    int pivotIndex = getPivot(indices, start, end);
    pivotIndex = partition(indices, start, end, pivotIndex);
    if (n == pivotIndex) return n;
    else if (n < pivotIndex) end = pivotIndex;
    else start = pivotIndex + 1;
  }
}

int KDTreeNode::partition(std::vector<int> &indices, int start, int end,
                          int pivot)
{
  double plane = points(indices[pivot], splitDim);
  swap(indices, pivot, end - 1);
  int i = start, j = end - 2;
  while (i < j) {
    if (points(indices[i], splitDim) >= plane) {
      swap(indices, i, j);
      j--;
    }
    else i++;
  }
  if (points(indices[i], splitDim) < plane) {
    swap(indices, i + 1, end - 1);
    return i + 1;
  }
  else {
    swap(indices, i, end - 1);
    return i;
  }
}

KDTreeNode::KDTreeNode(const Eigen::MatrixXd &points, std::vector<int> &indices,
                       int start, int end, int dim, int maxDim) : points(points)
{
  if (end - start == 1) {
    index = indices[start];
    splitDim = -1;
    leftChild = NULL;
    rightChild = NULL;
    return;
  }

  splitDim = dim;
  do {
    int middle = findMedian(indices, start, end);
    splitPlane = points(indices[middle], splitDim);
    while (middle > 0 && points(indices[middle - 1], splitDim) == splitPlane) {
      middle--;
    }
    for (int i = start; i < middle; i++) {
      if (points(indices[i], splitDim) >= splitPlane) std::cout << "Trouble left " << points(indices[i], splitDim) << " " << splitPlane << " " << splitDim << " " << start << " " << end << " " << " " << middle << " " << i << "\n";
    }
    for (int i = middle; i < end; i++) {
      if (points(indices[i], splitDim) < splitPlane) std::cout << "Trouble right\n";
    }

    if (middle > start) {
      leftChild = new KDTreeNode(points, indices, start, middle,
                                 (dim + 1) % maxDim, maxDim);
      rightChild = new KDTreeNode(points, indices, middle, end,
                                 (dim + 1) % maxDim, maxDim);
      return;
    }
    splitDim = (splitDim + 1) % maxDim;
  } while (splitDim != dim);
  //can only get here if all points in range are equal, just make it a leaf
  index = indices[start];
  splitDim = -1;
  leftChild = NULL;
  rightChild = NULL;
}

KDTreeNode::~KDTreeNode()
{
  if (leftChild != NULL) delete leftChild;
  if (rightChild != NULL) delete rightChild;
}

int KDTreeNode::lookup(Eigen::VectorXd &query)
{
  if (splitDim == -1) { //leaf node
    double bestDist = pow(points(index, 0), 2) +
                      pow(points(index, 1), 2) +
                      pow(points(index, 2), 2);

    return index;
  }
  if (leftChild == NULL) return rightChild->lookup(query);
  int closestInd = query[splitDim] < splitPlane ?
                   leftChild->lookup(query) :
                   rightChild->lookup(query);

  double bestDist = pow(points(closestInd, 0), 2) +
                    pow(points(closestInd, 1), 2) +
                    pow(points(closestInd, 2), 2);

  if (pow(query[splitDim] - splitPlane, 2) < bestDist) {
    int otherInd = query[splitDim] < splitPlane ?
                   leftChild->lookup(query) :
                   rightChild->lookup(query);
    double otherDist = pow(points(otherInd, 0), 2) +
                       pow(points(otherInd, 1), 2) +
                       pow(points(otherInd, 2), 2);
    if (otherDist < bestDist) return otherInd;
    else return closestInd;
  }
  return closestInd;
}

KDTree::KDTree(const Eigen::MatrixXd &points)
{
  k = points.cols();
  std::vector<int> indices(points.rows());
  for (int i = 0; i < points.rows(); i++) indices[i] = i;
  root = new KDTreeNode(points, indices, 0, points.rows(), 0, k);
}

int KDTree::lookup(Eigen::VectorXd &query)
{
  return root->lookup(query);
}

KDTree::~KDTree()
{
  delete root;
}
