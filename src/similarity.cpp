#include <Rcpp.h>
using namespace Rcpp;
#include <algorithm>    // std::set_intersection, std::sort
#include <vector>       // std::vector

//Get set intersection
std::vector<int> intersect(std::vector<int> s1, std::vector<int> s2){
  // Sizes of both the sets
  int size_s1 = s1.size();
  int size_s2 = s2.size();
  //v holds the intersection set
  //allocate it as big as both sets combined just in case
  std::vector<int> v(size_s1 + size_s2);
  //iterator for use in std::set_intersection()
  std::vector<int>::iterator it;

  //sort the input vectors
  std::sort(s1.begin(), s1.end());
  std::sort(s2.begin(), s2.end());

  // Get the intersection set
  std::vector<int>::iterator itr = std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), v.begin());
  //v will have the intersection and then zeros
  //delete the zeros
  v.resize(itr-v.begin());

  return v;
}


//Calculate size of set intersection
int size_intersect(std::vector<int> s1, std::vector<int> s2){

  //get intersection set
 std::vector<int> v = intersect(s1, s2);
  return v.size();
}

// Calculate Jaccard similarity of two sets
// [[Rcpp::export]]
double jaccard_index(NumericVector s1, NumericVector s2){
  int size_s1 = s1.size();
  int size_s2 = s2.size();
  std::vector<int> s1_std = Rcpp::as<std::vector<int>>(s1);
  std::vector<int> s2_std = Rcpp::as<std::vector<int>>(s2);
  int size_in = size_intersect(s1_std, s2_std);
  // Calculate the Jaccard index
  // size of intersection/size of union
  double jaccard_in = size_in / (size_s1 + size_s2 - size_in);

  // Return the Jaccard index
  return jaccard_in;
}



//Function to enumerate children of a given node
std::vector<int> get_children_std(int this_node, std::vector<int> tree_nodes, std::vector<int> tree_parents){
  std::vector<int> children; //holds immediate children of this_node -- starting length 0
  //count number of nodes whose parent is this_node, if any
  int n_children = std::count(tree_parents.begin(), tree_parents.end(), this_node);
  if(n_children>0){
    //iterate over all nodes and *find* those whose parent is this_node
    for (int i=0; i<tree_nodes.size(); i++){
      if(tree_parents[i] == this_node){
        //add to children
        children.push_back(tree_nodes[i]);
      }
    }
  }

  return children; //note that children may be length-0 if there are no children
}

//Rcpp-format wrapper for get_children_std
// [[Rcpp::export]]
IntegerVector get_children(int this_node, IntegerVector tree_nodes, IntegerVector tree_parents){
  std::vector<int> nodes_std = Rcpp::as<std::vector<int>>(tree_nodes);
  std::vector<int> parents_std = Rcpp::as<std::vector<int>>(tree_parents);
  std::vector<int> children_std = get_children_std(this_node, nodes_std, parents_std);
  //convert to IntegerVector before return
  return Rcpp::wrap(children_std);
}

//Function to enumerate all descendants of each of a given list of nodes
std::vector<std::vector<int>> get_descendants_std(std::vector<int> these_nodes, std::vector<int> tree_nodes, std::vector<int> tree_parents){
  int n = these_nodes.size(); //length of these_nodes
  std::vector<std::vector<int>> descendants(n); //n elements, each itself an int vector of initial size zero
//iterate over these_nodes and get their children
for(int i=0; i<n; i++){
 std::vector<int> children = get_children_std(these_nodes[i], tree_nodes, tree_parents);
  //children will be size 0 if we have reached a tip. Otherwise:
  if(children.size() > 0){
    //append the children to the i-th descendants vector
  descendants[i].insert(descendants[i].end(), children.begin(), children.end());
    //now call get_descendants only on *one* child at a time
    for(int j=0; j<children.size(); j++){
      //make a 1-element vector with this child to pass to get_descendants
    std::vector<int> this_child{children[j]};
      //recursively call get_descendants on this_child
    std::vector<std::vector<int>> new_desc = get_descendants_std(this_child, tree_nodes, tree_parents);
    //new_desc will be a 1-element vector of an int vector
    //extract the int vector and append it to descendants[i]
    descendants[i].insert(descendants[i].end(), new_desc[0].begin(), new_desc[0].end() );
    } //end for j-0, j<n, j++
  } //end if children.size()>0
} //end for i=0, i<n, i++
return descendants;
}

//Rcpp-format wrapper for get_descendants_std
// [[Rcpp::export]]
List get_descendants(IntegerVector these_nodes, IntegerVector tree_nodes, IntegerVector tree_parents){
  std::vector<int> these_nodes_std = Rcpp::as<std::vector<int>>(these_nodes);
  std::vector<int> nodes_std = Rcpp::as<std::vector<int>>(tree_nodes);
  std::vector<int> parents_std = Rcpp::as<std::vector<int>>(tree_parents);
  std::vector<std::vector<int>> desc_std = get_descendants_std(these_nodes_std, nodes_std, parents_std);
  //convert to List before return
  return Rcpp::wrap(desc_std);
}

//Function to calculate information content for a set of nodes, with respect to a given tree
//Here, information content is a function of number of descendants of each node, relative to total tree size
std::vector<double> calc_IC_std(std::vector<int> these_nodes,
                                std::vector<int> tree_nodes,
                                std::vector<int> tree_parents,
                                NumericMatrix information_content){
  int n = these_nodes.size();
  //int tree_size = tree_nodes.size();

  //get descendants for each node
  //std::vector<std::vector<int>> descendants(n);
  //descendants = get_descendants_std(these_nodes, tree_nodes, tree_parents);

  std::vector<double> IC(n);

  //For each node: count descendants & calculate IC accordingly
  //for (int i=0; i<n; i++){
  //  int n_desc = descendants[i].size();
  //  IC[i] = 1 - log(1 + n_desc)/log(tree_size);
  //}

  for (int i = 0; i < n; i++){
    for (int j = 0; j < information_content.nrow(); j++){
      if (these_nodes[i] == information_content(j, 0)){
        IC[i] = information_content(j, 4);
        break;
      }
    }
  }

  return IC;
}

//Function to calculate information content of nodes with reference to a hierarchical taxonomy
// [[Rcpp::export]]
NumericVector calc_IC(IntegerVector these_nodes, //node IDs for which to calc IC
                      IntegerVector tree_nodes, //node IDs for whole taxonomy
                      IntegerVector tree_parents, //parent IDs for each node
                      NumericMatrix information_content){
  std::vector<int> these_nodes_std = Rcpp::as<std::vector<int>>(these_nodes);
  std::vector<int> nodes_std = Rcpp::as<std::vector<int>>(tree_nodes);
  std::vector<int> parents_std = Rcpp::as<std::vector<int>>(tree_parents);
  int n = these_nodes.size();
  std::vector<double> IC(n);
  IC = calc_IC_std(these_nodes_std, nodes_std, parents_std, information_content);
  return Rcpp::wrap(IC);
}



//
//Function to get parent of a given node
std::vector<std::vector<int>> get_parent_std(std::vector<int> these_nodes, std::vector<int> tree_nodes, std::vector<int> tree_parents){
  int n = these_nodes.size();   //number of nodes to check
  std::vector<std::vector<int>> parent(n); //declare a "list" to hold the parent of each node. Each list element will be either size 0 or size 1
  for(int i=0; i<n; i++){
    for(int j = 0; j<tree_nodes.size(); j++){
      if(tree_nodes[j] == these_nodes[i]){
        //check if this parent is in tree_nodes; if not, it's root and has no parent
        //(this works because the root node will be listed with some placeholder parent like 0, -1, or even NA
        if(std::count(tree_nodes.begin(), tree_nodes.end(), tree_parents[j])>0){
          //handle the potential case when root node is listed with itself as parent:
          //this won't happen if nodes and parents are generated with get_tree_df(),
          //but might happen in another situation
          if(tree_parents[j] != tree_nodes[j]){
            //add this parent to the list of parents
            parent[i].push_back(tree_parents[j]);
          }
        }
      }
    }
  }
  return parent;
}

//Rcpp wrapper for get_parent_std
// [[Rcpp::export]]
List get_parents(IntegerVector these_nodes, IntegerVector tree_nodes, IntegerVector tree_parents){
  std::vector<int> these_nodes_std = Rcpp::as<std::vector<int>>(these_nodes);
  std::vector<int> nodes_std = Rcpp::as<std::vector<int>>(tree_nodes);
  std::vector<int> parents_std = Rcpp::as<std::vector<int>>(tree_parents);
  std::vector<std::vector<int>> parent_std = get_parent_std(these_nodes_std, nodes_std, parents_std);
  //convert to List before return
  return Rcpp::wrap(parent_std);
}

//
//Function to get all ancestors of a given node
std::vector<std::vector<int>> get_ancestors_std(std::vector<int> these_nodes, std::vector<int> tree_nodes, std::vector<int> tree_parents){
  int n = these_nodes.size();
  std::vector<std::vector<int>> ancestors(n);
  for(int i=0; i<n; i++){
    //pack this_node into a one-element vector
    std::vector<int> this_node{these_nodes[i]};
    //get a one-element "list" of its parent, if any
    std::vector<std::vector<int>> this_parent = get_parent_std(this_node, tree_nodes, tree_parents);
    //if there was any parent, i.e. if we have not reached root:
    if(this_parent[0].size() > 0){
    //extract the parent and append it to ancestors
    ancestors[i].insert(ancestors[i].end(), this_parent[0].begin(), this_parent[0].end());
     //call get_ancestors_std on this_parent
     std::vector<std::vector<int>> new_anc = get_ancestors_std(this_parent[0], tree_nodes, tree_parents);
    //new_anc will be a 1-element vector of an int vector
    //extract the int vector and insert it to ancestors[i]
        ancestors[i].insert(ancestors[i].end(), new_anc[0].begin(), new_anc[0].end() );
      }
  }

  return ancestors;
}

//Rcpp wrapper for get_ancestors_std
//Do this because Rcpp Vector classes aren't very efficient with regards to .insert or .push_back methods,
//so it's better to do the actual operations on std::vector instead
// [[Rcpp::export]]
List get_ancestors(IntegerVector these_nodes, IntegerVector tree_nodes, IntegerVector tree_parents){
  std::vector<int> these_nodes_std = Rcpp::as<std::vector<int>>(these_nodes);
  std::vector<int> nodes_std = Rcpp::as<std::vector<int>>(tree_nodes);
  std::vector<int> parents_std = Rcpp::as<std::vector<int>>(tree_parents);
  std::vector<std::vector<int>> ancestors_std = get_ancestors_std(these_nodes_std, nodes_std, parents_std);
  //convert to IntegerVector before return
  return Rcpp::wrap(ancestors_std);
}

//Function to get most recent common ancestor of two nodes wrt a hierarchical taxonomy
std::vector<int> get_MRCA_std(int node1, int node2, std::vector<int> tree_nodes, std::vector<int> tree_parents){
  //get ancesors of node_1 and node_2
  std::vector<int> node1_v{node1};
  std::vector<int> node2_v{node2};
  std::vector<std::vector<int>> anc1_v = get_ancestors_std(node1_v, tree_nodes, tree_parents);
  std::vector<std::vector<int>> anc2_v = get_ancestors_std(node2_v, tree_nodes, tree_parents);
  //create vectors of ancestors with the nodes themselves at the beginning
  //this means that if node1 is an ancestor of node2,
  //that their MRCA will be node1 itself (or vice versa if node2 is ancestor of node1).
  std::vector<int> anc1{node1};
  anc1.insert(anc1.end(), anc1_v[0].begin(), anc1_v[0].end());
  std::vector<int> anc2{node2};
  anc2.insert(anc2.end(), anc2_v[0].begin(), anc2_v[0].end());


  //iterator pointing to the first element in anc1 that matches any element in anc2
  std::vector<int>::iterator itr = std::find_first_of(anc1.begin(),
                                                      anc1.end(),
                                                      anc2.begin(),
                                                      anc2.end());

  std::vector<int> MRCA;
  if(itr!=anc1.end()){ //if there was any match at all
    MRCA.push_back(*itr); //append it to the MRCA output vector
  }
  //otherwise if no match, MRCA output vector remains empty.
  //this *should* never happen if we have a rooted tree, but who knows.

  return MRCA;


}

//Rcpp wrapper for get_MRCA_std
//Do this because Rcpp Vector classes aren't very efficient with regards to .insert or .push_back methods,
//so it's better to do the actual operations on std::vector instead
// [[Rcpp::export]]
IntegerVector get_MRCA(int node1, int node2, IntegerVector tree_nodes, IntegerVector tree_parents){
  std::vector<int> nodes_std = Rcpp::as<std::vector<int>>(tree_nodes);
  std::vector<int> parents_std = Rcpp::as<std::vector<int>>(tree_parents);
  std::vector<int> MRCA = get_MRCA_std(node1, node2, nodes_std, parents_std);
  return Rcpp::wrap(MRCA);
}

//Function to calculate the Resnik similarity of two nodes
double get_resnik_std(int node1, int node2, std::vector<int> tree_nodes, std::vector<int> tree_parents,
                      NumericMatrix information_content){
  //this is IC of the MRCA
  std::vector<int> MRCA(1);
  MRCA = get_MRCA_std(node1, node2, tree_nodes, tree_parents);
  std::vector<double> resnik(1);
  for (int i = 0; i < information_content.nrow(); i++){
    if (information_content(i,0) == MRCA[0]){
      resnik[0] = information_content(i, 4);
      break;
    }
  }
  //resnik = calc_IC_std(MRCA, tree_nodes, tree_parents);
  return resnik[0]; //resnik will only have 1 element, for the MRCA, so return it
}

//Rcpp wrapper for get_resnik_std
//Do this because Rcpp Vector classes aren't very efficient with regards to .insert or .push_back methods,
//so it's better to do the actual operations on std::vector instead
// [[Rcpp::export]]
double get_resnik(int node1, int node2, IntegerVector tree_nodes, IntegerVector tree_parents, NumericMatrix information_content){
  std::vector<int> tree_nodes_std = Rcpp::as<std::vector<int>>(tree_nodes);
  std::vector<int> tree_parents_std = Rcpp::as<std::vector<int>>(tree_parents);
  double resnik = get_resnik_std(node1, node2, tree_nodes_std, tree_parents_std, information_content);
  return resnik;
}

//Function to calculate the Lin similarity of two nodes
double get_lin_std(int node1, int node2, std::vector<int> tree_nodes, std::vector<int> tree_parents, NumericMatrix information_content){
  double resnik = get_resnik_std(node1, node2, tree_nodes, tree_parents, information_content);
  std::vector<double> IC(2);
  std::vector<int>node12{node1, node2};
  IC = calc_IC_std(node12, tree_nodes, tree_parents, information_content);
  double lin = (2*resnik)/(IC[0] + IC[1]);
  return lin;
}

//Rcpp wrapper for get_lin_std
//Do this because Rcpp Vector classes aren't very efficient with regards to .insert or .push_back methods,
//so it's better to do the actual operations on std::vector instead
// [[Rcpp::export]]
double get_lin(int node1, int node2, IntegerVector tree_nodes, IntegerVector tree_parents, NumericMatrix information_content){
  std::vector<int> tree_nodes_std = Rcpp::as<std::vector<int>>(tree_nodes);
  std::vector<int> tree_parents_std = Rcpp::as<std::vector<int>>(tree_parents);
  double lin = get_lin_std(node1, node2, tree_nodes_std, tree_parents_std, information_content);
  return lin;
}

//Function to calculate the Jiang & Conrath similiarty of two nodes
double get_jiang_conrath_std(int node1, int node2, std::vector<int> tree_nodes, std::vector<int> tree_parents, NumericMatrix information_content){
  double resnik = get_resnik_std(node1, node2, tree_nodes, tree_parents, information_content);
  std::vector<double> IC(2);
  std::vector<int>node12{node1, node2};
  IC = calc_IC_std(node12, tree_nodes, tree_parents, information_content);
  double jiang_conrath = 1 - (IC[0] + IC[1] - 2*resnik)/2;

  return jiang_conrath;
}

//Rcpp wrapper for get_jiang_conrath_std
//Do this because Rcpp Vector classes aren't very efficient with regards to .insert or .push_back methods,
//so it's better to do the actual operations on std::vector instead
// [[Rcpp::export]]
double get_jiang_conrath(int node1, int node2, IntegerVector tree_nodes, IntegerVector tree_parents, NumericMatrix information_content){
  std::vector<int> tree_nodes_std = Rcpp::as<std::vector<int>>(tree_nodes);
  std::vector<int> tree_parents_std = Rcpp::as<std::vector<int>>(tree_parents);
  double jiang_conrath = get_jiang_conrath_std(node1, node2, tree_nodes_std, tree_parents_std, information_content);
  return jiang_conrath;
}

//Function to calculate the Jaccard index of two nodes
double get_jaccard_std(int node1, int node2, std::vector<int> tree_nodes, std::vector<int> tree_parents){
  //put each node into a 1-element vector
  std::vector<int> node_v{node1, node2};
  //get ancestors of each node
  std::vector<std::vector<int>> anc_v = get_ancestors_std(node_v, tree_nodes, tree_parents);
  //calc size of intersection of ancestors
  int size_in = size_intersect(anc_v[0], anc_v[1]);
  // Calculate the Jaccard index
  // size of intersection/size of union
  double jaccard = size_in / (anc_v[0].size() + anc_v[1].size() - size_in);
  return jaccard;
}

//Rcpp wrapper for get_jaccard_std
//Do this because Rcpp Vector classes aren't very efficient with regards to .insert or .push_back methods,
//so it's better to do the actual operations on std::vector instead
// [[Rcpp::export]]
double get_jaccard(int node1, int node2, IntegerVector tree_nodes, IntegerVector tree_parents){
  std::vector<int> tree_nodes_std = Rcpp::as<std::vector<int>>(tree_nodes);
  std::vector<int> tree_parents_std = Rcpp::as<std::vector<int>>(tree_parents);
  double jaccard = get_jaccard_std(node1, node2, tree_nodes_std, tree_parents_std);
  return jaccard;
}

// Function to return the pairwise similarities of two lists of nodes
// [[Rcpp::export]]
NumericMatrix get_similarity(IntegerVector nodes1, //node IDs in set 1
                             IntegerVector nodes2, //node IDs in set 2
                             IntegerVector tree_nodes, //all nodes of tree
                             IntegerVector tree_parents, //all parents of tree
                             int sim_metric, // which similarity metric: 1 = Jaccard, 2 = Resnik, 3 = Lin, 4 = Jiang and Conrath
                             NumericMatrix information_content,
                             bool upper_tri = true){

  int n1 = nodes1.size();
  int n2 = nodes2.size();
  //std::vector versions of nodes
  //Do this because Rcpp Vector classes aren't very efficient with regards to .insert or .push_back methods,
  //so it's better to do the actual operations on std::vector instead
  std::vector<int> nodes1_std = Rcpp::as<std::vector<int>>(nodes1);
  std::vector<int> nodes2_std = Rcpp::as<std::vector<int>>(nodes2);
  std::vector<int> tree_nodes_std = Rcpp::as<std::vector<int>>(tree_nodes);
  std::vector<int> tree_parents_std = Rcpp::as<std::vector<int>>(tree_parents);
  //declare similarity matrix to store outputs
  NumericMatrix m(n1, n2);

  if (upper_tri){
  for(int i = 0; i < n1; ++i){
    for(int j = i; j < n2; ++j){
      int node_i = nodes1_std[i];
      int node_j = nodes2_std[j];
      switch(sim_metric){
      case 1: //Jaccard similarity
        m(i,j) = get_jaccard_std(node_i, node_j, tree_nodes_std, tree_parents_std);
        break;
      case 2: //Resnik similarity
        m(i,j) = get_resnik_std(node_i, node_j, tree_nodes_std, tree_parents_std, information_content);
        break;
      case 3: //Lin similarity
        m(i,j) = get_lin_std(node_i, node_j, tree_nodes_std, tree_parents_std, information_content);
        break;
      case 4: //Jiang and Conrath similarity
        m(i,j) = get_jiang_conrath_std(node_i, node_j, tree_nodes_std, tree_parents_std, information_content);
        break;
      }
    }
  }
  }
  else
  {
    for(int i = 0; i < n1; ++i){
      for(int j = 0; j < n2; ++j){
        int node_i = nodes1_std[i];
        int node_j = nodes2_std[j];
        switch(sim_metric){
        case 1: //Jaccard similarity
          m(i,j) = get_jaccard_std(node_i, node_j, tree_nodes_std, tree_parents_std);
          break;
        case 2: //Resnik similarity
          m(i,j) = get_resnik_std(node_i, node_j, tree_nodes_std, tree_parents_std, information_content);
          break;
        case 3: //Lin similarity
          m(i,j) = get_lin_std(node_i, node_j, tree_nodes_std, tree_parents_std, information_content);
          break;
        case 4: //Jiang and Conrath similarity
          m(i,j) = get_jiang_conrath_std(node_i, node_j, tree_nodes_std, tree_parents_std, information_content);
          break;
        }
      }
    }
  }
  return m;
}
