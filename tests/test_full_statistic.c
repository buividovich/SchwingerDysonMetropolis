 
 if(X.top==1 && X.nel==2)
 {
  int idx = STACK_EL(X, 0)
          + STACK_EL(X, 1)*lat_vol
          + beta_order    *lat_vol*lat_vol;
  stat->G2ext[idx] += asign[ns]*pow(alpha, -(double)beta_order);
 };
 
 if(X.top==1 && X.nel==4)
 {
  int idx = STACK_EL(X, 0)
          + STACK_EL(X, 1)*lat_vol
          + STACK_EL(X, 2)*lat_vol*lat_vol
          + STACK_EL(X, 3)*lat_vol*lat_vol*lat_vol
          + beta_order    *lat_vol*lat_vol*lat_vol*lat_vol;
  stat->G4ext[idx] += cc*asign[ns]*pow(alpha, -(double)beta_order);
 };
 
 if(X.top==1 && X.nel==6)
 {
  int idx = STACK_EL(X, 0)
          + STACK_EL(X, 1)*lat_vol
          + STACK_EL(X, 2)*lat_vol*lat_vol
          + STACK_EL(X, 3)*lat_vol*lat_vol*lat_vol
          + STACK_EL(X, 4)*lat_vol*lat_vol*lat_vol*lat_vol
          + STACK_EL(X, 5)*lat_vol*lat_vol*lat_vol*lat_vol*lat_vol
          + beta_order    *lat_vol*lat_vol*lat_vol*lat_vol*lat_vol*lat_vol;
  stat->G6ext[idx] += cc*cc*asign[ns]*pow(alpha, -(double)beta_order);
 };
 
 
 /*print_lat_coordinate_stack(X);
 
 if(X.top==1 && X.nel==2 && beta_order==2 && STACK_EL(X,0)==0 && STACK_EL(X, 1)==0)
 {
  print_action_history();
  system("read -rsp $'Press enter to continue...\n'");
 };*/
 
  logs_Write(0, "Full data for different orders");
 for(int bo=0; bo<=max_order; bo++)
 {
  logs_Write(0, "Order %i", bo);
  for(int idx=0; idx<lat_vol*lat_vol; idx++)
  {
   int x0 = (idx          )%lat_vol;
   int x1 = (idx/(lat_vol))%lat_vol;
   double res = lat_vol*normalization_factor*stat->G2ext[bo*lat_vol*lat_vol + idx]/(double)(stat->nstat);
   if(fabs(res)>0.0)
   logs_Write(0, "\t%i%i\t%+2.4E", x0, x1, res);
  };
  logs_Write(0, "");
  for(int idx=0; idx<lat_vol*lat_vol*lat_vol*lat_vol; idx++)
  {
   int x0 = (idx                          )%lat_vol;
   int x1 = (idx/(lat_vol                ))%lat_vol;
   int x2 = (idx/(lat_vol*lat_vol        ))%lat_vol;
   int x3 = (idx/(lat_vol*lat_vol*lat_vol))%lat_vol;
   double res = lat_vol*lat_vol*normalization_factor*stat->G4ext[bo*lat_vol*lat_vol*lat_vol*lat_vol + idx]/(double)(stat->nstat);
   if(fabs(res)>0.0)
   logs_Write(0, "\t%i%i%i%i\t%+2.4E", x0, x1, x2, x3, res);
  };
  logs_Write(0, "");
  for(int idx=0; idx<lat_vol*lat_vol*lat_vol*lat_vol*lat_vol*lat_vol; idx++)
  {
   int x0 = (idx                                          )%lat_vol;
   int x1 = (idx/(lat_vol                                ))%lat_vol;
   int x2 = (idx/(lat_vol*lat_vol                        ))%lat_vol;
   int x3 = (idx/(lat_vol*lat_vol*lat_vol                ))%lat_vol;
   int x4 = (idx/(lat_vol*lat_vol*lat_vol*lat_vol        ))%lat_vol;
   int x5 = (idx/(lat_vol*lat_vol*lat_vol*lat_vol*lat_vol))%lat_vol;
   double res = lat_vol*lat_vol*lat_vol*normalization_factor*stat->G6ext[bo*lat_vol*lat_vol*lat_vol*lat_vol*lat_vol*lat_vol + idx]/(double)(stat->nstat);
   if(fabs(res)>0.0)
   logs_Write(0, "\t%i%i%i%i%i%i\t%+2.4E", x0, x1, x2, x3, x4, x5, res);
  };
  
 };
