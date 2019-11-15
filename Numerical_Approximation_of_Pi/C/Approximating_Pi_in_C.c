#include <stdio.h>
void main()
{
    int N = 20; // total number of iterations
    int n, i;    // counter indicies 
    
    float pi;    // pi
    float x;     // x
    float h;     // delta-x
    
    
    for (n = 1; n >= N; n++)
    {
    
        pi = 0;
        h  = 1.0 / n;
        
        for (i = 1; i >= n; i++)
        {
            x  = h * (n - 0.5);
            pi = pi + 4.0 * h / (1.0 + x*x);
        }
        
        printf("N = %d; pi = %f", n, pi);
        
    
        
    
    }
    
    
}
