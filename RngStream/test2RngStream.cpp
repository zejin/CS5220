/*  Programme pour tester le generateur   RngStream.cpp         */

#include <iostream>
#include "RngStream.h"
using namespace std;

int main ()
{
   double sum = 0.0;
   int  i;
   RngStream g1 ("g1");
   RngStream g2 ("g2");
   RngStream g3 ("g3");

   cout.precision(13);
   cout << "Initial states of g1, g2, and g3:\n\n";
   g1.WriteState ();   g2.WriteState ();   g3.WriteState ();
   sum = g2.RandU01 () + g3.RandU01 ();
   for (i = 0;  i < 12345; i++)
      g2.RandU01 ();

   g1.AdvanceState (5, 3);   
   cout << "State of g1 after advancing by 2^5 + 3 = 35 steps:\n";
   g1.WriteState ();
   cout << "g1.RandU01 () = " << g1.RandU01 () << "\n\n";
 
   g1.ResetStartStream ();
   for (i = 0;  i < 35; i++)    g1.AdvanceState (0,1);
   cout << "State of g1 after reset and advancing 35 times by 1:\n";
   g1.WriteState ();
   cout << "g1.RandU01 () = " << g1.RandU01 () << "\n\n";
 
   g1.ResetStartStream ();
   long sumi = 0;
   for (i = 0;  i < 35; i++)    sumi += g1.RandInt (1, 10);
   cout << "State of g1 after reset and 35 calls to RandInt (1, 10):\n";
   g1.WriteState ();
   cout << "   sum of 35 integers in [1, 10] = " << sumi << "\n\n";
   sum += sumi / 100.0;
   cout << "g1.RandU01 () = " << g1.RandU01 () << "\n\n";

   double sum3 = 0.0;
   g1.ResetStartStream ();
   g1.IncreasedPrecis (true);
   sumi = 0;
   for (i = 0;  i < 17; i++)     sumi += g1.RandInt (1, 10);
   cout << "State of g1 after reset, IncreasedPrecis (true) and 17 calls"
        << " to RandInt (1, 10):\n";
   g1.WriteState ();
   g1.IncreasedPrecis (false);
   g1.RandInt (1, 10);
   cout << "State of g1 after IncreasedPrecis (false) and 1 call to RandInt\n";
   g1.WriteState ();
   sum3 = sumi / 10.0;

   g1.ResetStartStream ();
   g1.IncreasedPrecis (true);
   for (i = 0;  i < 17; i++)    sum3 += g1.RandU01 ();
   cout << "State of g1 after reset, IncreasedPrecis (true) and 17 calls"
        <<   " to RandU01:\n";
   g1.WriteState ();
   g1.IncreasedPrecis (false);
   g1.RandU01 ();
   cout << "State of g1 after IncreasedPrecis (false) and 1 call to RandU01\n";
   g1.WriteState ();
   sum += sum3 / 10.0;

   sum3 = 0.0;
   cout << "Sum of first 100 output values from stream g3:\n";
   for (i = 0;  i < 100;  i++) {
      sum3 += g3.RandU01 ();
   }
   cout << "   sum = " << sum3 << "\n\n";
   sum += sum3 / 10.0;

   cout << "\nReset stream g3 to its initial seed.\n";
   g3.ResetStartStream ();
   cout << "First 5 output values from stream g3:\n";
   for (i=1; i<=5; i++)
      cout << g3.RandU01 () << "\n";
   sum += g3.RandU01 ();

   cout << "\nReset stream g3 to the next SubStream, 4 times.\n";
   for (i=1; i<=4; i++)
      g3.ResetNextSubstream ();
   cout << "First 5 output values from stream g3, fourth SubStream:\n";
   for (i=1; i<=5; i++)
      cout << g3.RandU01 () << "\n";
   sum += g3.RandU01 ();

   cout << "\nReset stream g2 to the beginning of SubStream.\n";
   g2.ResetStartSubstream ();
   cout << " Sum of 100000 values from stream g2 with double precision:   ";
   sum3 = 0.0;
   g2.IncreasedPrecis (true);
   for (i=1; i<=100000; i++)
      sum3 += g2.RandU01 ();
   cout << sum3 << "\n";
   sum += sum3 / 10000.0;
   g2.IncreasedPrecis (false);

   g3.SetAntithetic (true);
   cout << " Sum of 100000 antithetic output values from stream g3:   ";
   sum3 = 0.0;
   for (i=1; i<=100000; i++)
      sum3 += g3.RandU01 ();
   cout << sum3 << "\n";
   sum += sum3 / 10000.0;

   cout << "\nSetPackageSeed to seed = { 1, 1, 1, 1, 1, 1 }\n";
   unsigned long germe[6] = { 1, 1, 1, 1, 1, 1 };
   RngStream::SetPackageSeed (germe);

   cout << "\nDeclare an array of 4 named streams"
        <<  " and write their full state\n";
   RngStream gar[4] = { "Poisson", "Laplace", "Galois", "Cantor" };
   for  (i = 0; i < 4; i++)
      gar[i].WriteStateFull ();

   cout << "Jump stream Galois by 2^127 steps backward\n" ;
   gar[2].AdvanceState (-127, 0);
   gar[2].WriteState ();
   gar[2].ResetNextSubstream ();

   for  (i = 0; i < 4; i++)
      sum += gar[i].RandU01 ();

   cout << "--------------------------------------\n";
   cout << "Final Sum = " << sum << "\n\n";

   return 0;
}
