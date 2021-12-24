#include <iostream>
#include <cstdlib>
#include <cmath>
#include <queue>


using namespace std;

void newGraph(double* A, int n, double* DA)//A - ishodnaya matrica, n - razmer, DA - vihodnoi graph
{
    queue<int> DX, DY;

    int cap = n*n;
    double front[cap];
    for(int i = 0; i < cap; i++)
    {
        DA[i] = 0; //5
        front[i] = INFINITY; //6
    }

    queue<int> Q1, Q2, Xfree, Yfree;//7

    double Xdouble[n], Ydouble[n];//dlya so4etanij
    for (int i = 0; i < n; i++)
    {
        Xdouble[i] = NAN;
        Ydouble[i] = NAN;
    }

    int x = 0, y = 0; 
    for(;x<n; x++)//8 pereproverit'
    {
        if (Xdouble[x] == NAN)//9
        {
            Q2.push(x);
            Xfree.push(x);//11
            DX.push(x);
            front[x]=0;//12
        }
    }

    while(!Q2.empty())
    {
        Q1.push(Q2.front);
        Q2.pop(); //15
    }

    int z = 0;
    while(!Q1.empty() || !Yfree.empty() || !Q2.empty())//16 + 37 opyat' ne uveren 4to tuda kuda nado zasunul
    {
        x = q1.front();//18
        for(; y <= n; y++)//19 pereproverit'
        {
            if(A[x + y*n] != 0 && front[x] < front[y])//20 
            {
                DA[x + y*n] = 1;//23 pust' budet poka bez vesa
                if (front[y] == INFINITY)//24
                {
                    DY.push(y);
                    z = Ydouble[y];//26
                    front[y] = front[x] + 1;//27
                    if (z!=0)//28 to4no 0? mojet vse taki nan?
                    {
                        DA[z + y*n] = 1;
                        DX.push(z);//30
                        front[z] = front[y] + 1;
                        Q2.push(z);//31
                    }
                    else Yfree.push(y);//33
                }
            }
        }
    }
}

void oldGraph(double* A, int n, double* DA) // legacy, for reference
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