#include <Rcpp.h>
using namespace Rcpp;
#include <algorithm>    // std::set_intersection, std::sort
#include <vector>       // std::vector

// Function to return the Jaccard index of two sets
// [[Rcpp::export]]
double jaccard_index(NumericVector s1, NumericVector s2)
{
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

  // Calculate the Jaccard index
  // size of intersection/size of union
  double jaccard_in = size_in
    / (size_s1 + size_s2 - size_in);

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
   //and list2
   for(int j = 0; j < length2; ++j){
      jaccard(i,j) = jaccard_index(list1[i], list2[j]);
   }
 }

 return jaccard;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
x <- 1:4
y <- 3:6
jaccard_index(x, y)
*/
