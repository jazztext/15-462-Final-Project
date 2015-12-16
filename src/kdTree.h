#include <vector>
#include <eigen3/Eigen/Dense>

class KDTreeNode
{
  public:
    const Eigen::MatrixXd &points;
    int splitDim;
    int index; //only for leaf nodes
    double splitPlane;
    KDTreeNode *leftChild;
    KDTreeNode *rightChild;

    KDTreeNode(const Eigen::MatrixXd &points, std::vector<int> &indices,
               int start, int end, int dim, int maxDim);
    int lookup(Eigen::VectorXd &query);
    ~KDTreeNode();

    bool operator()(int i, int j);

  private:
    int partition5(std::vector<int> &indices, int start, int end);
    int getPivot(std::vector<int> &indices, int start, int end);
    int findMedian(std::vector<int> &indices, int start, int end);
    int partition(std::vector<int> &indices, int start, int end, int pivot);


};

class KDTree
{
  public:
    //takes points to be inserted as nxk matrix
    KDTree(const Eigen::MatrixXd &points);
    int lookup(Eigen::VectorXd &query);

    ~KDTree();

  protected:
    int k;
    KDTreeNode *root;
};
