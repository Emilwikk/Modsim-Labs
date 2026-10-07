#include <stdio.h>
#include <math.h>

int main()
{
  double c;

  printf("Please provide a floating point value! ");
  scanf("%lf", &c);

  printf("square root: %g\n", sqrt(c));

  return 0;
}
