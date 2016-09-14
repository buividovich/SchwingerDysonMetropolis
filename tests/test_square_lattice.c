
void test_are_neighbors()
{
 logs_Write(0, "TESTING THE ARE_NEIGHBORS ROUTINE");
 char S1[64], S2[64], SD[32];
 for(int x1=0; x1<lat_vol; x1++)
 {
  int nc = 0;
  for(int x2=0; x2<lat_vol; x2++)
  {
   int dir = are_neighbors(x1, x2);
   if(dir>=0)
   {
    sprintf_coords(S1, x1);
    sprintf_coords(S2, x2);
    sprintf(SD, "%s", (dir<DIM? "+" : "-"));
    if(dir>=DIM)
     dir -= DIM;
    logs_Write(0, "\t %s<->%s are neighbors, direction is %s%i", S1, S2, SD, dir);
    nc ++;
   };
  };
  printf("%i neighbors in total\n", nc);
 }; 
}; 

/* logs_Write(0, "TESTING THE LAT_DISTANCE ROUTINE");
 char S1[64], S2[64], SD[32];
 for(int x1=0; x1<lat_vol; x1++)
  for(int x2=0; x2<lat_vol; x2++)
  {
   int xd = lat_distance(x1, x2);
   sprintf_coords(S1, x1);
   sprintf_coords(S2, x2);
   sprintf_coords(SD, xd);
   logs_Write(0, "\t %s<->%s\t %s", S1, S2, SD);
  };
 
 return EXIT_SUCCESS;*/
