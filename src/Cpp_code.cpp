//// Cpp functions 2015-3-11 Thong Pham
#include <Rcpp.h>
#include <iostream>    
#include <algorithm>    
#include <vector>      

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

double my_zeroin(double ax,double bx,std::function <double (double)> f,double tol,long max_iter)  ;

// [[Rcpp::export(".normalized_constant")]]
int normalized_constant(      NumericVector& norm, 
                        const NumericMatrix& degree, 
                        const NumericVector& theta,
                        const NumericVector& f, 
                        const NumericMatrix& offset_tk,
                        const double&        offset
                        ) {
    long T = degree.nrow();     // number of time-steps
    long N = degree.ncol();     // number of nodes
    long K = offset_tk.ncol();  // maximum degree
    #pragma omp parallel for
    for (long i = 0; i < T; i++) {
        double total = 0; 
        for (long j = 0; j < N; j++)
            if (degree(i,j) >= 0) {
                total += theta(degree(i,j)) * f(j);
            }
        for (long k = 0; k < K; ++k) {
            //printf("%f",offset_tk(i,k));
            total += offset * offset_tk(i,k) * theta(k);
        }
        norm(i) = total;
    }
    return 0;
}

// [[Rcpp::export(".normalized_constant_alpha")]]
int normalized_constant_alpha(      NumericVector& norm,
                              const double       & alpha,   
                              const double       & PA_offset,  
                              const NumericMatrix& degree, 
                              const NumericVector& theta,
                              const NumericVector& f, 
                              const NumericMatrix& offset_tk,
                              const double&        offset) {
  long T = degree.nrow();     // number of time-steps
  long N = degree.ncol();     // number of nodes
  long K = offset_tk.ncol();  // maximum degree 
  //printf("alpha is: %f \n",alpha);
  //#pragma omp parallel for
  for (long i = 0; i < T; i++) {
    double total = 0;
    for (long j = 0; j < N; j++)
      if (degree(i,j) >= 0) {
          if  (degree(i,j) > 0)  {
              total += (pow(theta.at(degree.at(i,j)),alpha)) * f.at(j);
          //printf("%f ",pow(theta.at(degree.at(i,j)),alpha));
          }
          else 
              total += PA_offset * f.at(j);  
    }
      //printf("%f \n",total);
    for (long k = 1; k < K; ++k)
        total += offset_tk.at(i,k) * pow(theta.at(k),alpha);
      //printf("%f \n",total);
    total += offset_tk.at(i,0) * PA_offset;
    norm.at(i) = total;
      //printf("%f \n",total);
      //printf("------------------------\n");
  }
  return 0;
}


// [[Rcpp::export(".get_stats")]]
int get_stats(CharacterVector    & time_stamp,
              CharacterVector    & unique_stamp,
              const NumericVector& in_node, 
              const NumericVector& out_node,
              const NumericVector& all_node,
              const NumericVector& ok_node,
              const NumericVector& bin_vector,
              const long max_node_id,
              const int  undirected,
              const int  only_PA,
              CharacterVector& time_vector,
              NumericVector& Sum_m_k, 
              NumericMatrix& n_tk,
              NumericVector& m_tk,
              NumericVector& m_t,
              NumericMatrix& offset_tk, 
              NumericVector& z_j, 
              NumericMatrix& node_degree,
              NumericMatrix& offset_m_tk,
              const int      only_true_deg,
              const long     deg_max,
              NumericVector& center_bin,
              NumericVector& appear_time) {
  long N     = all_node.size(); 
  long N_new = ok_node.size();
  long T          = time_vector.size();
  long K          = n_tk.ncol();
  //std::cout << "Begin the program /n";
  
  std::vector<long>      node_array(max_node_id + 1,0); //the index used when indexing all arrays whose length is N
  std::vector<int>       ok_array(max_node_id + 1,0);   //boolean check whether a node is used or not  
  std::vector<int>       is_appear(N,0);
  std::vector<int>       appear_onestep(N,0);
  std::vector<long>      ok_index(max_node_id + 1,0);   //the index used when indexing all arrays whose length is N_new
  std::vector<long>      degree_vector(N,-1);
  std::vector<long>      degree_vector_onestep(N,-1);
  std::vector<long>      n_tk_vector(K,0);
  std::vector<long>      m_tk_vector(K,0);
  std::vector<double>    count_bin(K,0);        // Number of degrees actually inside that bin
  //std::vector<double>  center_bin(K,0);       // Logarithmic center of the bin
  std::vector<long>      is_in_bin(deg_max,0);  // Mark a degree as appeared in bin
  std::vector<long>      z_j_vector(N_new,0);
  std::vector<long>      offset_tk_vector(K,0);
  std::vector<long>      offset_m_tk_vector(K,0);
  
  //std::cout << "start first loop here /n";
  for (long long i = 0; i < N; ++ i) {
      node_array[all_node(i)] = i; 
  }
  //std::cout << "finished first loop here /n";
  
  for (long long i = 0; i < N_new; ++ i) {
      ok_array[ok_node(i)] = 1; 
      ok_index[ok_node(i)] = i;
     
  }
  //std::cout << "Passed here /n";
  
  long t  = 0;
  long edge_count = 0;

  while ((t < T)) {
      checkUserInterrupt();  
      if (t > 0)  {
          if (0 == only_PA) { 
              for (long long j = 0; j < ok_node.size(); ++ j) {
                  if (degree_vector.at(node_array.at(ok_node(j))) >= 0) {
                      node_degree.at(t - 1,ok_index.at(ok_node(j))) = bin_vector(degree_vector.at(node_array.at(ok_node(j))));
                  }else {
                      node_degree.at(t - 1,ok_index.at(ok_node(j))) = -1;
                  }
              }
          }
          if (0 == only_true_deg)
              for (long k = 0; k < K; ++k) {
                  n_tk(t - 1,k)          = n_tk_vector.at(k);
                  if (0 == only_PA) {  
                      offset_tk.at(t - 1,k)   = offset_tk_vector.at(k);
                      offset_m_tk.at(t - 1,k) = offset_m_tk_vector.at(k);  
                      offset_m_tk_vector.at(k) = 0;
                  }
              }
      }    
      if (0 == only_true_deg)
          for (long k = 0; k < K; ++k)
              m_tk_vector.at(k) = 0;
      
      if (0 == only_true_deg)
          if (0 == only_PA)      
              for (long long j = 0; j < N_new; ++j)
                  z_j_vector.at(j) = 0;    
      
      while ((edge_count < in_node.size()) && (time_stamp(edge_count) == time_vector(t))) {
           long in_node_id  = in_node(edge_count);
           long out_node_id = out_node(edge_count);
           //std::cout << out_node_id << " " << in_node_id << " / ";
           // if the node is an isolated node
           // Assumption: an isolated node only appears one time
           if (-1 == in_node_id || -1 == out_node_id) {
             //std::cout << "in isolated node when crash / ";   
              if (-1 != in_node_id) {
                  long in_node_ind = node_array.at(in_node(edge_count));   
                
                  if (0 == appear_onestep.at(in_node_ind)) {    
                      degree_vector.at(in_node_ind)  = 0;
                      //degree_vector_onestep.at(in_node_ind) = 0;
                  }
                  
                  if (0 == only_true_deg)   
                      if (0 == is_appear[in_node_ind])  
                          if ((0 == only_PA) && (1 == ok_array.at(in_node(edge_count))))  
                              appear_time.at(ok_index.at(in_node(edge_count))) = t + 1; 
                  
                  if (0 == only_true_deg)
                      if (0 == appear_onestep.at(in_node_ind))         
                          ++n_tk_vector.at(bin_vector(degree_vector.at(in_node_ind))); 
                      
                  if (0 == only_true_deg)
                      if (0 == appear_onestep.at(in_node_ind))         
                          if ((0 == only_PA) && (0 == ok_array.at(in_node(edge_count))))   
                              ++offset_tk_vector.at(bin_vector(degree_vector.at(in_node_ind)));   
                      
                  //is_appear[in_node_ind]         = 1;  
                  appear_onestep.at(in_node_ind) = 1;
              }
              if (-1 != out_node_id) {
                  long out_node_ind = node_array.at(out_node(edge_count));  
                  
                  if (0 == appear_onestep.at(out_node_ind)) {
                      degree_vector.at(out_node_ind)  = 0;
                      //degree_vector_onestep.at(out_node_ind) = 0;
                  }
                  
                  if (0 == only_true_deg)   
                      if ((0 == only_PA) && (1 == ok_array.at(out_node_ind)))  
                          appear_time.at(ok_index.at(out_node_ind)) = t + 1; 
                      
                  if (0 == only_true_deg)
                      if (0 ==  appear_onestep.at(out_node_ind))  
                          ++n_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));
                      
                  if (0 == only_true_deg)
                      if (0 ==  appear_onestep.at(out_node_ind))  
                          if ((0 == only_PA) && (0 == ok_array.at(out_node_ind)))   
                              ++offset_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));    
                      
                  //is_appear[out_node_ind]         = 1;  
                  appear_onestep.at(out_node_ind) = 1;    
              }
           } else {
             //std::cout << "in normal node when crash / ";
             //std::cout << out_node(edge_count) << " ";
             //std::cout << in_node(edge_count) << " / ";
             
             
             //std::cout << "Before converted ID: ";
             long in_node_ind  = node_array.at(in_node(edge_count)); 
             long out_node_ind = node_array.at(out_node(edge_count));   
               
           
          
             //std::cout << out_node(edge_count) << " ";
             //std::cout << in_node(edge_count) << " / ";
             
             //std::cout << " Before convert: ";
             
             //std::cout << "Appear of in node: " << is_appear[in_node_ind] << " " << appear_onestep[in_node_ind] << ". ";
             
             //std::cout << "Appear of out node: " << is_appear[out_node_ind] << " " << appear_onestep[out_node_ind] << ". ";
          
           // consider the in-node first
           //the in-node has not appeared in the previous time step   
           if (0 == is_appear[in_node_ind]) {
                 if (0 == only_true_deg)   
                     if ((0 == only_PA) && (1 == ok_array.at(in_node(edge_count))))  
                         appear_time.at(ok_index.at(in_node(edge_count))) = t + 1;  
              //the in-node has not appeared in previous edges of the current time step
              if (0 == appear_onestep.at(in_node_ind)) { 
                  appear_onestep.at(in_node_ind)  = 1;
                  degree_vector.at(in_node_ind)   = 1;
                  //appear_onestep_in.at(in_node_ind)  = 1;
                  if (0 == only_true_deg)
                      ++n_tk_vector.at(bin_vector(degree_vector.at(in_node_ind))); 
                  if (0 == only_true_deg)
                      if ((0 == only_PA) && (0 == ok_array.at(in_node(edge_count))))   
                          ++offset_tk_vector.at(bin_vector(degree_vector.at(in_node_ind))); 
              
              }
              else { //the in-node has already appeared in some edges of the current time step
                  if (0 == only_true_deg)    
                      if ((0 == only_PA)&& (0 == ok_array.at(in_node(edge_count))))    
                          --offset_tk_vector.at(bin_vector(degree_vector.at(in_node_ind))); 
                  if (0 == only_true_deg)
                      --n_tk_vector.at(bin_vector(degree_vector.at(in_node_ind))); 
                  ++degree_vector.at(in_node_ind); 
                  if (0 == only_true_deg)
                      ++n_tk_vector.at(bin_vector(degree_vector.at(in_node_ind))); 
                  if (0 == only_true_deg)
                      if ((0 == only_PA) && (0 == ok_array.at(in_node(edge_count))))   
                          ++offset_tk_vector.at(bin_vector(degree_vector.at(in_node_ind))); 
              }     
        } // the in-node is already appeared in the previous time-step
        else { 
            //std::cout << "inside of in_node loop for already appeared node. ";      
          
            if (0 == only_true_deg)  
                if ((0 == only_PA) && (1 == ok_array.at(in_node(edge_count))))
                    z_j(ok_index.at(in_node(edge_count)))++;
                
            //std::cout << "reach 1 ";   
            if (0 == only_true_deg)    {
                //std::cout << "in_node_ind: " << in_node_ind;
                //std::cout << ". degree one_step: " << degree_vector_onestep.at(in_node_ind);
                //std::cout << " . degree: " << degree_vector.at(in_node_ind);
                //std::cout << ". bin_vector: " << bin_vector(degree_vector_onestep.at(in_node_ind)) ;
                //std::cout << ". m_tk: " << m_tk_vector.at(bin_vector(degree_vector_onestep.at(in_node_ind))) << ". "; 
                ++m_tk_vector.at(bin_vector(degree_vector_onestep.at(in_node_ind)));
                 //std::cout << "reach 1.1 ";   
                  if (0 == is_in_bin.at(degree_vector_onestep.at(in_node_ind))) {
                   // std::cout << "reach 1.2 ";   
                      is_in_bin.at(degree_vector_onestep.at(in_node_ind)) = 1;
                     // std::cout << "reach 1.3 ";     
                      if (degree_vector_onestep.at(in_node_ind) > 0) {
                       // std::cout << "reach 1.4 ";   
                          center_bin.at(bin_vector(degree_vector_onestep.at(in_node_ind))) += log10((double) degree_vector_onestep.at(in_node_ind));
                         // std::cout << "reach 1.5 ";   
                          count_bin.at(bin_vector(degree_vector_onestep.at(in_node_ind)))  += 1;
                          //std::cout << "reach 1.6 ";   
                      }
                  }
            }
            // std::cout << "reach 2 ";
            
            if (0 == only_true_deg)
               if ((0 == only_PA) && (0 == ok_array.at(in_node(edge_count)))) {
                   ++offset_m_tk_vector.at(bin_vector(degree_vector_onestep.at(in_node_ind)));    
                   if (0 == is_in_bin.at(degree_vector_onestep.at(in_node_ind))) {
                       is_in_bin.at(degree_vector_onestep.at(in_node_ind)) = 1;
                       if (degree_vector_onestep.at(in_node_ind) > 0) {
                           center_bin.at(bin_vector(degree_vector_onestep.at(in_node_ind))) += log10( (double) degree_vector_onestep.at(in_node_ind));
                           count_bin.at(bin_vector(degree_vector_onestep.at(in_node_ind)))  += 1;
                       }
                 }
               }
            //   std::cout << "reach 3 ";      
            if (0 == only_true_deg)
                --n_tk_vector.at(bin_vector(degree_vector.at(in_node_ind)));
            
            // std::cout << "reach 4 ";
            
            if (0 == only_true_deg)
                if ((0 == only_PA) && (0 == ok_array.at(in_node(edge_count)))) {                  
                    -- offset_tk_vector.at(bin_vector(degree_vector.at(in_node_ind)));        
               }
                
            // std::cout << "reach 5 ";       
           ++degree_vector.at(in_node_ind);
           
           // std::cout << "reach 6 ";
           
           if (0 == only_true_deg)      
               ++n_tk_vector.at(bin_vector(degree_vector.at(in_node_ind))); 
           
           // std::cout << "reach 7 ";
           
           if (0 == only_true_deg)
               if ((0 == only_PA) && (0 == ok_array.at(in_node(edge_count)))) {                       
                   ++offset_tk_vector.at(bin_vector(degree_vector.at(in_node_ind)));    
               }
        }
      //consider next the Out node
      //the out node has not appeared in the previous time step
      if (0 == is_appear.at(out_node_ind)) {
           // std::cout << "outside of out_node loop. ";    
           if (0 == only_true_deg)   
              if ((0 == only_PA) && (1 == ok_array.at(out_node(edge_count))))  
                   appear_time.at(ok_index.at(out_node(edge_count))) = t + 1;   
          //the network is undirected, so this out_node is also counted
          if (1 == undirected) {
             // std::cout << "reach inside here. ";  
              if (0 == appear_onestep.at(out_node_ind)) {
                  degree_vector.at(out_node_ind)   = 1;  
                   if (0 == only_true_deg) 
                       ++n_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));  
                   if (0 == only_true_deg)
                       if ((0 == only_PA) && (0 == ok_array.at(out_node(edge_count))))
                           ++offset_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));     
              }
              else {
                  if (0 == only_true_deg)  
                      if ((0 == only_PA) && (0 == ok_array.at(out_node(edge_count)))) 
                          --offset_tk_vector.at(bin_vector(degree_vector.at(out_node_ind))); 
                  if (0 == only_true_deg)    
                      --n_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));    
                  degree_vector.at(out_node_ind)   = 1;   
                  if (0 == only_true_deg)
                      ++n_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));   
                  if (0 == only_true_deg)
                      if ((0 == only_PA) && (0 == ok_array.at(out_node(edge_count)))) 
                          ++offset_tk_vector.at(bin_vector(degree_vector.at(out_node_ind))); 
              }
          }
          //the network is directed, so this out_node is not counted
          else {
               if (0 == appear_onestep.at(out_node_ind)) { 
                  degree_vector.at(out_node_ind)   = 0; 
                  if (0 == only_true_deg)
                      ++n_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));
                  if (0 == only_true_deg)
                      if ((0 == only_PA) && (0 == ok_array.at(out_node(edge_count)))) {                       
                          ++offset_tk_vector.at(bin_vector(degree_vector.at(out_node_ind))); 
                      }
              }    
          }
          appear_onestep.at(out_node_ind)       = 1;
      } 
      // the out node is already appeared in the previous time-step
      else {
           //the network is undirected, so this out_node is also counted  
          if (1 == undirected) {  
              if (0 == only_true_deg)  
                  if ((0 == only_PA) && (0 == ok_array.at(out_node(edge_count)))) {                  
                      -- offset_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));        
                  }
              if (0 == only_true_deg)    
                  if ((0 == only_PA) && (1 == ok_array.at(out_node(edge_count))))
                      z_j.at(ok_index.at(out_node(edge_count)))++;
              if (0 == only_true_deg) {    
                  ++m_tk_vector.at(bin_vector(degree_vector_onestep.at(out_node_ind))); 
                  if (0 == is_in_bin.at(degree_vector_onestep.at(out_node_ind))) {
                      is_in_bin.at(degree_vector_onestep.at(out_node_ind)) = 1;
                      if (degree_vector_onestep.at(out_node_ind) > 0) {
                          center_bin.at(bin_vector(degree_vector_onestep.at(out_node_ind))) += log10((double) degree_vector_onestep.at(out_node_ind));
                          count_bin.at(bin_vector(degree_vector_onestep.at(out_node_ind)))  += 1;
                      }
                }
              }
              if (0 == only_true_deg)
                  if ((0 == only_PA) && (0 == ok_array.at(out_node(edge_count)))) {
                      ++offset_m_tk_vector.at(bin_vector(degree_vector_onestep.at(out_node_ind)));
                      if (0 == is_in_bin.at(degree_vector_onestep.at(out_node_ind))) {
                          is_in_bin.at(degree_vector_onestep.at(out_node_ind)) = 1;
                          if (degree_vector_onestep.at(out_node_ind) > 0) {
                              center_bin.at(bin_vector(degree_vector_onestep.at(out_node_ind))) += log10((double) degree_vector_onestep.at(out_node_ind));
                              count_bin.at(bin_vector(degree_vector_onestep.at(out_node_ind)))  += 1;
                      }
                    }
                  }
              if (0 == only_true_deg)    
                  --n_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));
              
              ++degree_vector.at(out_node_ind);
              
              if (0 == only_true_deg)
                  ++n_tk_vector.at(bin_vector(degree_vector.at(out_node_ind))); 
              
              if (0 == only_true_deg)
                  if ((0 == only_PA) && (0 == ok_array.at(out_node(edge_count)))) {                  
                      -- offset_tk_vector.at(bin_vector(degree_vector.at(out_node_ind)));        
                  }
          }
      }
      }
      ++edge_count; 
     }
  if (0 == only_true_deg)
      if (t > 0)  {
          for (long k = 0; k < K; ++k) {
              m_tk(t - 1,k)     = m_tk_vector.at(k);
              Sum_m_k(k)       += m_tk_vector.at(k);
              m_t(t - 1)       += m_tk_vector.at(k);
          
          }
          if (0 == only_PA) { 
              for (long i = 0; i < N_new; ++i) {
                  z_j(i) += z_j_vector.at(i);
              } 
          for (long k = 0; k < K; ++k) {
              offset_m_tk(t - 1,k) = offset_m_tk_vector.at(k);      
          }
          }
      }
     t++;
    
     for (long n = 0; n < (long) degree_vector.size(); ++n)
         degree_vector_onestep.at(n) = degree_vector.at(n);
     for (long i = 0; i < (long) is_appear.size(); ++i) {
         is_appear.at(i) = appear_onestep.at(i); 
         //appear_onestep.at(i)     = 0;
         //appear_onestep_in.at(i)  = 0;
         //appear_onestep_out.at(i) = 0;
     }
}
   for (unsigned long k = 0; k < count_bin.size(); ++ k)
       if (count_bin.at(k) != 0)  { 
           center_bin.at(k) /= count_bin.at(k);  
           center_bin.at(k)  = pow(10,center_bin.at(k));
       } 
       else center_bin.at(k) = 0; 
  //std::cout << "Outside loop " ;    
   return 0;
}
// [[Rcpp::export(".update_f")]]
int update_f(      NumericVector& f, 
             const NumericVector& non_zero_f,
             const NumericMatrix& degree, 
             const NumericVector& theta, 
             const NumericVector& z_j,
             const NumericVector& normalized_const, 
             const NumericVector& m_t, 
             const double         shape, 
             const double         rate,
             const double         offset) {
    long T        = degree.nrow();        // number of time-steps
    long N_nozero = non_zero_f.size();   // number of nodes
    #pragma omp parallel for
    for (long j = 0; j < N_nozero; j++) {
        double total = 0;
        for (long i = 0; i < T; i++)
            if ((degree(i,non_zero_f(j) - 1) >= 0) && (normalized_const(i) != 0)) {
                total += m_t(i) / normalized_const(i) * theta(degree(i,non_zero_f(j) - 1));
            }
    if (z_j(non_zero_f(j) - 1) + shape - 1 <= 0)
        f(non_zero_f(j) - 1) = f(non_zero_f(j) - 1);
    else 
        f(non_zero_f(j) - 1) = (z_j(non_zero_f(j) - 1) + shape - 1) / (total + rate);
    }
    return 0;
}



// [[Rcpp::export(".update_offset")]]
double update_offset(
                   const NumericMatrix& offset_n_tk,
                   const NumericMatrix& offset_m_tk, 
                   const NumericVector& theta, 
                   const NumericVector& normalized_const, 
                   const NumericVector& m_t, 
                   const double         shape, 
                   const double         rate) {
  long T        = offset_n_tk.nrow();        // number of time-steps
  long K        = offset_n_tk.ncol();
  double total1 = 0;
  double total2 = 0;
  double offset = 1;  
  #pragma omp parallel for reduction(+:total1, total2) 
  for (long i = 0; i < T; i++) {
      for (long k = 0; k < K; k++) {  
      if (normalized_const(i) != 0) 
          total1 += m_t(i) / normalized_const(i) * offset_n_tk(i,k) * theta(k);
      total2 += offset_m_tk(i,k);  
      }
  }
  //std::cout << "total 2: " << total2 << "; Total 1: \n" << total1 << "\n";
  //std::cout << "shape: " << shape << "; rate: \n" << rate << "\n";
  if (total2 + shape - 1 > 0)
      offset = (total2 + shape - 1.0)/(total1 + rate);
  //printf("%f ",offset);
  
  //std::cout << "offset in C: " << offset << "\n";
  return offset;
}

// // [[Rcpp::export(".update_offset_alpha")]]
// double update_offset_alpha( 
//                          const double       & alpha,       
//                          const NumericMatrix& offset_n_tk,
//                          const NumericMatrix& offset_m_tk, 
//                          const NumericVector& theta, 
//                          const NumericVector& normalized_const, 
//                          const NumericVector& m_t, 
//                          const double         shape, 
//                          const double         rate) {
//   long T        = offset_n_tk.nrow();        // number of time-steps
//   long K        = offset_n_tk.ncol();
//   double total1 = 0;
//   double total2 = 0;
//   double offset = 0;
// #pragma omp parallel for reduction(+:total1, total2) 
//   for (long i = 0; i < T; i++) {
//     for (long k = 0; k < K; k++)  
//       if (normalized_const(i) != 0) {
//         total1 += m_t(i) / normalized_const(i) * offset_n_tk(i,k) * pow(theta.at(k),alpha);
//         total2 += offset_m_tk(i,k);  
//       }
//   }
//   offset = (total2 + shape - 1)/(total1 + rate);
//   
//   return offset;
// }


// [[Rcpp::export(".update_f_alpha")]]
int update_f_alpha(      NumericVector& f, 
                   const NumericVector& non_zero_f,
                   const double       & alpha,
                   const double       & PA_offset,
                   const NumericMatrix& degree, 
                   const NumericVector& theta, 
                   const NumericVector& z_j,
                   const NumericVector& normalized_const, 
                   const NumericVector& m_t, 
                   const double         shape, 
                   const double         rate
) {
  long T        = degree.nrow();        // number of time-steps
  long N_nozero = non_zero_f.size();    // number of nodes
  //long N_nozero = degree.ncol();    // number of nodes
  //  printf("Alpha inside: %f \n",alpha);
   #pragma omp parallel for
   for (long j = 0; j < N_nozero; j++) {
     double total = 0;
     for (long i = 0; i < T; i++)
         if ((degree.at(i,non_zero_f(j) - 1) >= 0) && (normalized_const.at(i) != 0)) {
             if (degree.at(i,j) > 0)
                 total += m_t.at(i) / normalized_const.at(i) * (pow(theta.at(degree.at(i,non_zero_f.at(j) - 1)),alpha)) ;
              else
                  total += m_t.at(i) / normalized_const.at(i) * (PA_offset);
         }
     if (z_j.at(non_zero_f(j) - 1) + shape - 1 <= 0)
         f.at(non_zero_f.at(j) - 1) = f.at(non_zero_f.at(j) - 1);
     else
         f.at(non_zero_f.at(j) - 1) = (z_j.at(non_zero_f.at(j) - 1) + shape - 1)/(total + rate);
   }
  return 0;
}


// [[Rcpp::export(".update_f_new")]]
int update_f_new(      NumericVector& f, 
                       const NumericVector& non_zero_f,
                       const NumericMatrix& degree, 
                       const NumericVector& theta, 
                       const NumericVector& z_j,
                       const NumericVector& normalized_const, 
                       const NumericVector& m_t, 
                       const double         shape, 
                       const double         rate,
                       const double         offset,
                       const NumericVector& weight_f) {
  long T        = degree.nrow();        // number of time-steps
  long N_nozero = non_zero_f.size();   // number of nodes
  #pragma omp parallel for
  for (long j = 0; j < N_nozero; j++) {
    double total = 0;
    for (long i = 0; i < T; i++)
        if ((degree(i,non_zero_f(j) - 1) >= 0) && (normalized_const(i) != 0)) {
            total += m_t(i) / normalized_const(i) * theta(degree(i,non_zero_f(j) - 1));
        }
    if (z_j(non_zero_f(j) - 1) + shape / weight_f.at(non_zero_f(j) - 1) - 1 <= 0)
        f(non_zero_f(j) - 1) = f(non_zero_f(j) - 1);
    else 
        f(non_zero_f(j) - 1) = (z_j(non_zero_f(j) - 1) + shape / weight_f.at(non_zero_f(j) - 1) - 1) / 
         (total + rate / weight_f.at(non_zero_f(j) - 1));
  }
  return 0;
}


// Using weighting of f
// [[Rcpp::export(".update_f_alpha_new")]]
int update_f_alpha_new(      NumericVector& f, 
                         const NumericVector& non_zero_f,
                         const double       & alpha,
                         const double       & PA_offset,
                         const NumericMatrix& degree, 
                         const NumericVector& theta, 
                         const NumericVector& z_j,
                         const NumericVector& normalized_const, 
                         const NumericVector& m_t, 
                         const double         shape, 
                         const double         rate,
                         const NumericVector& weight_f) {
  long T        = degree.nrow();        // number of time-steps
  long N_nozero = non_zero_f.size();    // number of nodes
  //long N_nozero = degree.ncol();    // number of nodes
  //  printf("Alpha inside: %f \n",alpha);
  #pragma omp parallel for
  for (long j = 0; j < N_nozero; j++) {
    double total = 0;
    for (long i = 0; i < T; i++)
      if ((degree.at(i,non_zero_f(j) - 1) >= 0) && (normalized_const.at(i) != 0)) {
        if (degree.at(i,non_zero_f(j) - 1) > 0)
          total += m_t.at(i) / normalized_const.at(i) * (pow(theta.at(degree.at(i,non_zero_f.at(j) - 1)),alpha)) ;
        else
          total += m_t.at(i) / normalized_const.at(i) * (PA_offset);
      }
    if (z_j.at(non_zero_f(j) - 1) + shape /weight_f.at(non_zero_f.at(j) - 1) - 1 <= 0)
        f.at(non_zero_f.at(j) - 1) = f.at(non_zero_f.at(j) - 1);
    else
        f.at(non_zero_f.at(j) - 1) = (z_j.at(non_zero_f.at(j) - 1) + 
                                        shape / weight_f.at(non_zero_f.at(j) - 1) - 1)/(total + 
                                        rate / weight_f.at(non_zero_f.at(j) - 1));
  }
  return 0;
}

// [[Rcpp::export(".update_alpha_fast")]]
double update_alpha_fast(
    const NumericVector& non_zero_theta,  
    const NumericVector& norm,
    const NumericVector& f, 
    const double       & PA_offset,
    const NumericVector& theta,
    const NumericMatrix& degree, 
    const NumericVector& m_t,
    const NumericVector& Sum_m_k,
    const NumericMatrix& offset_tk,
    const double&        offset,
    const double         alpha_old
) {
  long T = degree.nrow();     // number of time-steps
  long N  = degree.ncol();     // number of nodes
  //long N2 = f.size();
  
  // printf("%d",N);
  // printf("%d",N2);
  long K = offset_tk.ncol();  // maximum degree 
  long length_theta = theta.size();
  double first_temp = 0;
  std::vector<double> coeff_degree(length_theta,0);

  #pragma omp parallel for \
  default(shared)          \
  reduction(+: first_temp) 
  for(long k = 0; k < Sum_m_k.size(); ++k) {
    if (theta.at(k) > 0) {
      first_temp += Sum_m_k.at(k) * log(theta.at(k)); //the log is due to taking log10, so we have to use 10-base   
    }
  }
  
  // #pragma omp parallel for 
  // for (long k = 1; k < K; ++k) {
  //     for (long t = 0; t < T; t++) {
  //         for (long i = 0; i < N; ++i)
  //             if (degree(t,i) == k && (theta.at(degree(t,i)) > 0) && (norm.at(t) > 0))  {
  //             //if (theta.at(degree(t,i)) > 0) 
  //             //all the logs here are due to take derivative, so we have to use natural base  
  //             coeff_degree.at(k) += m_t.at(t) / norm.at(t) * f.at(i) * log(theta.at(k)); 
  //         }
  //         if (theta.at(k) > 0 && (norm.at(t) > 0)) {
  //             coeff_degree.at(k) += m_t.at(t) / norm.at(t) * offset_tk(t,k) * log(theta.at(k));
  //               //printf("%f ",coeff_degree.at(k));
  //         }    
  //     }
  // }
  
  for (long t = 0; t < T; t++) {
      for (long i = 0; i < N; ++i)
          if (degree(t,i) > 0 && (theta.at(degree(t,i)) > 0) && (norm.at(t) > 0))  {
              //if (theta.at(degree(t,i)) > 0) 
                  //all the logs here are due to take derivative, so we have to use natural base  
                  coeff_degree.at(degree(t,i)) += m_t.at(t) / norm.at(t) * f.at(i) * log(theta.at(degree(t,i))); 
      }
      for (long k = 1; k < K; ++k)
          if (theta.at(k) > 0 && (norm.at(t) > 0)) {
              coeff_degree.at(k) += m_t.at(t) / norm.at(t) * offset_tk(t,k) * log(theta.at(k));
              //printf("%f ",coeff_degree.at(k));
          }
  }
  
  //printf("\n %f ",first_temp);
  auto f_1 = [&](double x) {
    double temp = 0;
    for(long k = 1; k < Sum_m_k.size(); ++k) {
        if (theta.at(k) >= 0) {
            if (theta.at(k) > 0)  
                temp += coeff_degree.at(k) * pow(theta.at(k),x);
        }
    }
    return(first_temp - temp); };
  double alpha = my_zeroin(-2, 2, f_1 , DBL_EPSILON,500);    
  //printf("alpha inside C: %f\n",alpha);    
  return alpha;
}


// [[Rcpp::export(".var_alpha")]]
double var_alpha(
    const double         alpha, 
    const NumericVector& non_zero_theta,  
    const NumericVector& norm,
    const NumericVector& f, 
    const double       & PA_offset,
    const NumericVector& theta,
    const NumericMatrix& degree, 
    const NumericVector& m_t,
    const NumericVector& Sum_m_k,
    const NumericMatrix& offset_tk,
    const double&        offset
) {
  long T = degree.nrow();     // number of time-steps
  long N  = degree.ncol();     // number of nodes
  //long N2 = f.size();
  
  // printf("%d",N);
  // printf("%d",N2);
  
  //long K = offset_tk.ncol();  // maximum degree 
  //long length_theta = theta.size();
  double var_alpha = 0; 
  

  for (long t = 0; t < T; t++) {
      double norm            = 0;
      double norm_derivative = 0;
      double norm_second_derivative = 0;
      for (long i = 0; i < N; ++i)
          if (degree(t,i) >= 0 && (theta.at(degree(t,i)) > 0))  {
              norm  += f.at(i) * pow(theta.at(degree(t,i)) , alpha);
              norm_derivative += pow(theta.at(degree(t,i)) , alpha) * log10((double) theta.at(degree(t,i))) * f.at(i); 
              norm_second_derivative +=  pow(theta.at(degree(t,i)) , alpha) * log10((double) theta.at(degree(t,i))) * log10((double) theta.at(degree(t,i))) * 
                                         f.at(i);  
          }
      double upper = norm_second_derivative * norm - norm_derivative * norm_derivative;   
      double lower = norm * norm;
      var_alpha -= m_t.at(t) * upper / lower;
  }
      
  return - 1.0 / var_alpha;
}

// // [[Rcpp::export(".update_PA_offset")]]
// double update_PA_offset(
//     const NumericVector& norm,
//     const NumericVector& f, 
//     const NumericMatrix& degree, 
//     const NumericVector& m_t,
//     const NumericVector& Sum_m_k,
//     const NumericMatrix& offset_tk) {
//   long T = degree.nrow();     // number of time-steps
//   long N = degree.ncol();     // number of nodes
//  //long K = offset_tk.ncol();  // maximum degree 
//   double second_temp = 0;
//   
//   for (long t = 0; t < T; t++) 
//     for (long i = 0; i < N; ++i) {
//       if ((degree(t,i) == 0) && (norm(t) != 0)){
//         second_temp += m_t.at(t)*f.at(i)/norm.at(t);
//       }
//     }
//   for (long t = 0; t < T; t++) 
//       second_temp += offset_tk(t,0) * m_t.at(t) /norm.at(t);
//   
//   double result = Sum_m_k.at(0) * 1.0 / second_temp;
//   return result;
// }

// [[Rcpp::export(".coeff_theta")]]
NumericVector coeff_theta( const NumericMatrix& degree,  
                           const NumericVector& f,
                           const NumericVector& normalized_const,  
                           const NumericVector& m_t, 
                           const int            length_theta
                           ) {
    int nrow = degree.nrow();
    int ncol = degree.ncol();
    NumericVector total(length_theta);
    //#pragma omp parallel 
    //{
    NumericVector total_temp(length_theta);
    for (int j = 0; j < length_theta; ++j) {
        total_temp[j] = 0;
    } 
    //#pragma omp for
    for (int j = 0; j < ncol; j++) {
        for (int i = 0; i < nrow; i++)
           if ((degree(i,j) >= 0) && (normalized_const(i) != 0)) {
                total[degree(i, j)] += f(j)* m_t(i) / normalized_const(i) ;
            }
    }
    //#pragma omp critical
    //{
    //for (int j = 0; j < length_theta; ++j)
    //    total[j] += total_temp[j];
    //}
   // }
    return total;
}

// [[Rcpp::export(".coeff_var")]]
NumericVector coeff_var(const NumericMatrix& degree,  
                        const NumericVector& f,
                        const NumericVector& normalized_const,
                        const NumericVector& m_t,
                        const NumericMatrix& offset, 
                        const int            length_theta) {
    int nrow = degree.nrow();
    int ncol = degree.ncol();
    NumericMatrix temp(nrow,length_theta);

    NumericVector total(length_theta);

    for (int j = 0; j < ncol; j++) {
        for (int t = 0; t < nrow; t++)
           if (degree(t,j) >= 0) {
                temp(t,degree(t,j)) += f(j);
            }
    }
    #pragma omp parallel for
    for (int k = 0; k < length_theta; k++) {
        for (int t = 0; t < nrow; t++)
            if (normalized_const[t] != 0)
                total(k) += pow(temp(t,k) + offset(t,k),2)* m_t(t) / pow(normalized_const(t),2);
    }
    return total;
}

// [[Rcpp::export(".cal_var_f")]]
int cal_var_f(      NumericVector& cov_f, 
              const NumericVector& non_zero_f,
              const NumericMatrix& degree,
              const NumericVector& theta,
              const NumericVector& f,
              const NumericVector& z_j,
              const NumericVector& normalized_const,
              const NumericVector& m_t, 
              const double         shape) {
    int T    = degree.nrow();
    int N    = non_zero_f.size();
    #pragma omp parallel for
    for (int j = 0; j < N; j++) {
        double total = 0;
        for (int i = 0; i < T; i++)
            if ((degree(i,non_zero_f(j) - 1) >= 0) && (normalized_const(i) != 0)) {
                total += m_t(i) / pow(normalized_const(i),2) * pow(theta(degree(i,non_zero_f(j) - 1)),2);
            }
          cov_f(j) = 1/(z_j(non_zero_f(j) - 1)/pow(f(non_zero_f(j) - 1),2) + - total +
                     (shape - 1)*pow(f(non_zero_f(j) - 1),2));
    }
    return 0;
}

// Using weighting of s
// [[Rcpp::export(".cal_var_f_new")]]
int cal_var_f_new(        NumericVector& cov_f, 
                    const NumericVector& non_zero_f,
                    const NumericMatrix& degree,
                    const NumericVector& theta,
                    const NumericVector& f,
                    const NumericVector& z_j,
                    const NumericVector& normalized_const,
                    const NumericVector& m_t, 
                    const double         shape,
                    const NumericVector& weight_f) {
  int T    = degree.nrow();
  int N    = non_zero_f.size();
  #pragma omp parallel for
  for (int j = 0; j < N; j++) {
    double total = 0;
    for (int i = 0; i < T; i++)
      if ((degree(i,non_zero_f(j) - 1) >= 0) && (normalized_const(i) != 0)) {
        total += m_t(i) / pow(normalized_const(i),2) * pow(theta(degree(i,non_zero_f(j) - 1)),2);
      }
    cov_f(j) = 1/(z_j(non_zero_f(j) - 1)/pow(f(non_zero_f(j) - 1),2) + - total +
        (shape / weight_f.at(non_zero_f.at(j) - 1) - 1) * pow(f(non_zero_f(j) - 1),2));
  }
  return 0;
}

double my_zeroin(double ax,double bx,std::function <double (double)> f,double tol,long max_iter)    /* An estimate to the root	*/
{
  double a,b,c;				/* Abscissae, descr. see above	*/
  double fa;				/* f(a)				*/
  double fb;				/* f(b)				*/
  double fc;				/* f(c)				*/
  
  a = ax;  b = bx;  fa = f(a);  fb = f(b);
  c = a;   fc = fa;
  long count = 0;
  for(;;)		/* Main iteration loop	*/
  { count++;
    //printf("%ld ",count);
    if (count > max_iter)
      return b;
    double prev_step = b-a;		/* Distance from the last but one*/
    /* to the last approximation	*/
    double tol_act;			/* Actual tolerance		*/
    double p;      			/* Interpolation step is calcu- */
    double q;      			/* lated in the form p/q; divi- */
    /* sion operations is delayed   */
    /* until the last moment	*/
    double new_step;      		/* Step at this iteration       */
    
    if( fabs(fc) < fabs(fb) )
    {                         		/* Swap data for b to be the 	*/
    a = b;  b = c;  c = a;          /* best approximation		*/
    fa=fb;  fb=fc;  fc=fa;
    }
    tol_act = 2*LDBL_EPSILON*fabs(b) + tol/2;
    new_step = (c-b)/2;
    
    if( fabs(new_step) <= tol_act || fb == (double)0 )
      return b;				/* Acceptable approx. is found	*/
    
    /* Decide if the interpolation can be tried	*/
    if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
    && fabs(fa) > fabs(fb) )	/* and was in true direction,	*/
    {					/* Interpolatiom may be tried	*/
    double t1,cb,t2;
      cb = c-b;
      if( a==c )			/* If we have only two distinct	*/
      {				/* points linear interpolation 	*/
    t1 = fb/fa;			/* can only be applied		*/
    p = cb*t1;
    q = 1.0 - t1;
      }
      else				/* Quadric inverse interpolation*/
      {
        q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
        p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
        q = (q-1.0) * (t1-1.0) * (t2-1.0);
      }
      if( p>(double)0 )		/* p was calculated with the op-*/
    q = -q;			/* posite sign; make p positive	*/
    else				/* and assign possible minus to	*/
    p = -p;			/* q				*/
    
    if( p < (0.75*cb*q-fabs(tol_act*q)/2)	/* If b+p/q falls in [b,c]*/
    && p < fabs(prev_step*q/2) )	/* and isn't too large	*/
    new_step = p/q;			/* it is accepted	*/
    /* If p/q is too large then the	*/
    /* bissection procedure can 	*/
    /* reduce [b,c] range to more	*/
    /* extent			*/
    }
    
    if( fabs(new_step) < tol_act )	/* Adjust the step to be not less*/
    {
    if ( new_step > (double)0 ) {	/* than tolerance		*/
    new_step = tol_act;
    } else {
          new_step = -tol_act;
    }
    }
    a = b;  fa = fb;			/* Save the previous approx.	*/
    b += new_step;  fb = f(b);	/* Do step to a new approxim.	*/
    if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) )
    {                 			/* Adjust c for it to have a sign*/
    c = a;  fc = fa;                  /* opposite to that of b	*/
    }
  }
}
