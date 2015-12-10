/*  Program to test the random number streams file:    RngStream.c   */

#include <stdio.h>
#include "RngStream.h"

int main (void)
{
   double sum;
   double sum3;
   long sumi;
   int  i;
   RngStream g1, g2, g3;
   RngStream gar[4];
   unsigned long germe[6] = { 1, 1, 1, 1, 1, 1 };

   g1 = RngStream_CreateStream ("g1");
   g2 = RngStream_CreateStream ("g2");
   g3 = RngStream_CreateStream ("g3");

   sum = RngStream_RandU01 (g2) + RngStream_RandU01 (g3);

   RngStream_AdvanceState (g1, 5, 3);   
   sum += RngStream_RandU01 (g1);

   RngStream_ResetStartStream (g1);
   for (i = 0;  i < 35; i++)
      RngStream_AdvanceState (g1, 0, 1);
   sum += RngStream_RandU01 (g1);

   RngStream_ResetStartStream (g1);
   sumi = 0;
   for (i = 0;  i < 35; i++)
      sumi += RngStream_RandInt (g1, 1, 10);
   sum += sumi / 100.0;

   sum3 = 0.0;
   for (i = 0;  i < 100;  i++) {
      sum3 += RngStream_RandU01 (g3);
   }
   sum += sum3 / 10.0;

   RngStream_ResetStartStream (g3);
   for (i=1; i<=5; i++)
      sum += RngStream_RandU01 (g3);

   for (i=0; i<4; i++)
      RngStream_ResetNextSubstream (g3);
   for (i=0; i<5; i++)
      sum += RngStream_RandU01 (g3);

   RngStream_ResetStartSubstream (g3);
   for (i=0; i<5; i++)
      sum += RngStream_RandU01 (g3);

   RngStream_ResetNextSubstream (g2);
   sum3 = 0.0;
   for (i=1; i<=100000; i++)
      sum3 += RngStream_RandU01 (g2);
   sum += sum3 / 10000.0;

   RngStream_SetAntithetic (g3, 1);
   sum3 = 0.0;
   for (i=1; i<=100000; i++)
      sum3 += RngStream_RandU01 (g3);
   sum += sum3 / 10000.0;

   RngStream_SetPackageSeed (germe);
   gar[0] = RngStream_CreateStream ("Poisson");
   gar[1] = RngStream_CreateStream ("Laplace");
   gar[2] = RngStream_CreateStream ("Galois");
   gar[3] = RngStream_CreateStream ("Cantor");

   for  (i = 0; i < 4; i++)
      sum += RngStream_RandU01 (gar[i]);

   RngStream_AdvanceState (gar[2], -127, 0);
   sum += RngStream_RandU01 (gar[2]);

   RngStream_ResetNextSubstream (gar[2]);
   RngStream_IncreasedPrecis (gar[2], 1);
   sum3 = 0.0;
   for  (i = 0; i < 100000; i++)
      sum3 += RngStream_RandU01 (gar[2]);
   sum += sum3 / 10000.0;

   RngStream_SetAntithetic (gar[2], 1);
   sum3 = 0.0;
   for  (i = 0; i < 100000; i++)
      sum3 += RngStream_RandU01 (gar[2]);
   sum += sum3 / 10000.0;
   RngStream_SetAntithetic (gar[2], 0);

   RngStream_IncreasedPrecis (gar[2], 0);

   for  (i = 0; i < 4; i++)
      sum += RngStream_RandU01 (gar[i]);

   RngStream_DeleteStream (g1);
   RngStream_DeleteStream (g2);
   RngStream_DeleteStream (g3);

   printf ("-------------------------------------------\n");
   printf ("This program should print   39.697547445251\n");
   printf ("Actual program result =     %.12f\n\n", sum);

   return 0;
}
