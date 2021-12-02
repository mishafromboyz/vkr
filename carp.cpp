#include <iostream>
#include <cstdlib>
#include <cmath>


using namespace std;

void graph(double* A, int n, double* DA)
{
    int Dx[n], Dy[n];
    int counterForDx, counterForDy = 0;

    int cap = n*n;
    int front[cap];
    for(int i = 0; i<cap; i++)
        front[i] = INFINITY;
    
    int q1[n], q2[n], xFree[n], yFree[n];
    for (int i = 0; i<n; i++)
    {
        q1[i] = 0;
        q1[i] = 0;
        xFree[i] = 0;
        yFree[i] = 0;
    }

    double Xdouble[n], Ydouble[n];
    for (int i = 0; i<n; i++)
    {
        Xdouble[i] = 0;
        Ydouble[i] = 0;
    }

    for (int i = 0; i < n; i++)
    {
        for(int j = 0; j < cap; j+=n)
            if(A[i+j] != 0)
            {
                j = cap;
                Xdouble[i] = A[i+j]; //???
            }
    }

    for (int i = 0; i < cap; i+=n)
    {
        for(int j = 0; j < n; j++)
            if(A[i+j] != 0)
            {
                j = n;
                Ydouble[i] = A[i+j];
            }
    }

    int counterForq2 = 0, counterForXFree = 0;
    for(int i = 0; i<n; i++)
    {
        if (Xdouble[i] == 0)
        {
            q2[counterForq2] = i;
            xFree[counterForXFree] = i;
            counterForq2++;
            counterForXFree++;
            
            Dx[counterForDx] = i;
            counterForDx++;
            front[i] = 0; //12
        }
    }

    for(int i = 0; i<n; i++)
    {
        q1[i] = q2[i];
        q2[i] = 0; //15
    }


    int z=0;
    counterForq2 = 0;
    int counterForYFree = 0;
    for(int i = 0; i < n || q2[i] != 0; q1[i] != 0) //16 + uslovia iz 37. ne ponimau 4to takoe yfree != nil
    {
        DX[i] = q1[i]; //18 ??? 4e takoe x
        for(int j = 0; j<n; j++)
        {
            if(A[i+j] = 1 && front[i] < front[j]) //??????? 4e takoe y? a mojet nado odin front dlya x vtoroi dlya y?
                {
                    DA[i+j] = 1;
                    if (front[j] = INFINITY) //24
                    {
                        Dy[counterForDy] = j;
                        counterForDy++;
                        z = Ydouble[j];
                        front[j] = front[x] + 1;
                        if (z != 0)
                        {
                            DA[z+j] =1; //26 ?????????
                            DX[counterForDx] = z;
                            counterForDx++;//30

                            front[z] = front[j] + 1;
                            q2[counterForq2] = z;
                            counterForq2++;
                        }
                        else
                        {
                            yFree[counterForYFree] = j;
                            counterForYfree++;
                        } //34
                    } //35
                } //36

        }// ????? kuda det' 37 stroky k i ili k j? dumau 4to k i no gotov podumat' eshe

    }

}

//podumat' nad counterami