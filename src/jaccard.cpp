#include <Rcpp.h>
using namespace Rcpp;
#include <algorithm>    // std::set_intersection, std::sort
#include <vector>       // std::vector


double size_intersect(NumericVector s1, NumericVector s2){
  // Sizes of both the sets
  double size_s1 = s1.size();
  double size_s2 = s2.size();
  //v holds the intersection set
  //allocate it as big as both sets combined just in case
  std::vector<int> v(size_s1 + size_s2);
  //iterator for use in std::set_intersection()
  std::vector<int>::iterator it;

  //sort the input vectors
  s1 = s1.sort();
  s2 = s2.sort();

  // Get the intersection set
  it = std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), v.begin());
  //v will have the intersection and then zeros
  //delete the zeros
  v.resize(it-v.begin());

  // Size of the intersection set
  double size_in = v.size();

  return size_in;
}

// Function to return the Jaccard index of two sets
// [[Rcpp::export]]
double jaccard_index(NumericVector s1, NumericVector s2)
{
  int size_s1 = s1.size();
  int size_s2 = s2.size();
  int size_in = size_intersect(s1, s2);
  // Calculate the Jaccard index
  // size of intersection/size of union
  double jaccard_in = size_in / (size_s1 + size_s2 - size_in);

  // Return the Jaccard index
  return jaccard_in;
}

// Function to return the pairwise Jaccard indexes of two lists
// [[Rcpp::export]]
NumericVector get_jaccard(List list1, List list2){
  double length1 = list1.size();
  double length2 = list2.size();
  NumericMatrix jaccard(length1, length2);

 //iterate over elements of list1
 for(int i = 0; i < length1; ++i){
   //and then iterate over elements of list2
   for(int j = 0; j < length2; ++j){
      jaccard(i,j) = jaccard_index(list1[i], list2[j]);
   }
 }

 return jaccard;
}


// Function to return the pairwise similarities of two lists
// NumericMatrix get_similarity(List anc1, //ancestors of nodesin set 1 as list object
//                              List anc2, //ancestors of nodes in set 2 as list object
//                              NumericVector nodes1, //node IDs in set 1
//                              NumericVector nodes2, //node IDs in set 2
//                              NumericVector ic1, //information content of nodes in set 1
//                              NumericVector ic2, //information content of nodes in set 2
//                              int sim_metric){ // which similarity metric: 1 = Jaccard, 2 = Resnik, 3 = Lin, 4 = Jiang and Conrath
//   //sort the input vectors of nodes
//   NumericVector s1 = nodes1.sort();
//   NumericVector s2 = nodes2.sort();
//
//   //sorted union of nodes
//   NumericVector allnodes = get_union(s1, s2);
//
//   //declare similarity matrix to store outputs
//   NumericMatrix m(s1.size(), s1.size());
//
//   for(int i = 0; i < allnodes.size(); ++i){
//     for(int j = i; j < allnodes.size(); ++j){
//       int node_i = allnodes[i];
//       int node_j = allnodes[j];
//       if(std::find(s1.begin(), s1.end(), node_i) != s1.end()){
//
//       }
//
//     }
//   }
// switch(sim_metric){
// case 1:
//
// }
//
// }
