#include <iostream>
#include <math.h>
#include <map>
#include <string>
#include <vector>

using namespace std;


//find maximum value and indices   
vector<int> findmax (float mat[][10], int n){
    
    int p = 0;
    int q = 0;
    float aij; 
    float max =0.; 
    for (int i = 0; i < n; i += 1){
        
        for (int j = 0; j < n; j += 1){
            
            if (i != j){
                aij = abs(mat[i][j]);
                
                if(aij > max){
                    max = aij;
                    p = i;
                    q = j;
            }            

            }
        }
    }

    vector <int> a;
    a.push_back(p);
    a.push_back(q);
    
    return a;
    
    
}


int main(int argc,char* argv[]){
    int n = 10;
    float h = 1./n;
    float A=1.;
    
    
    float max = 0.;
    int l;
    int k;
    
    float t;
    float c;
    float s;
    
    float eps = pow(10.0, -4.0);
    

    float matrix[n][10] = {0};

//fills the matrix

    for (int i = 0; i<=(n-1); i +=1){
        matrix[i][i-1] = -1./pow(h,2);
        matrix[i][i] = 2./pow(h,2);
        matrix[i][i+1] = -1./pow(h,2);        
    }

    for (int i = 2; i<=(n-1); i+=1){
        matrix[0][i] = 0.;
    }
   
   

//Jacobi's method 

    while (A>eps){
        vector<int> a = findmax(matrix, n);

        l = a[0];
        k = a[1];
        
        float tau = (matrix[l][l]-matrix[k][k])/(2.*matrix[k][l]);
        
        if(tau >= 0){
            t = (-tau -sqrt(1+tau*tau));
        }else{
            t = (-tau +sqrt(1+tau*tau));   
        }
        
        c = 1./(sqrt(1+t*t));
        s = c*t;

        for (int i = 0; i < n; i+=1){
            if((i != k) && (i != l)){
                float ik = matrix[i][k];
                float il = matrix[i][l];
                matrix[i][k] = ik*c-il*s;
                matrix[i][l] = il*c+ik*s;
                matrix[k][i] = matrix[i][k];
                matrix[l][i] = matrix[i][l];
            }
        }

        float x = matrix[k][k]*pow(c,2.)-2.*matrix[k][l]*c*s+matrix[l][l]*pow(s, 2.);
        float y = matrix[l][l]*pow(c,2.)+2.*matrix[k][l]*c*s+matrix[k][k]*pow(s, 2.);

        matrix[k][l] = 0;
        matrix[l][k] = matrix[k][l];

        matrix[k][k] = x;
        matrix[l][l] = y;

        float z = 0;

        for(int u = 0; u < n; u++){
            for(int v = 0; v < n; v++){
                if(u != v){
                    z += pow(matrix[u][v],2);
                }

            }
        }
        
        A=sqrt(z);
       

      
        
    }

for (int i = 0; i < n; i++){
    cout << "Computed values: " << matrix[i][i] << endl; 
    }




//prints the matrix  
/*    
for (int i = 0; i<=(n-1); i++){
    for(int j = 0; j<=(n-1); j+=1){
        cout << matrix[i][j] << "  "; 
        }
    cout << endl;
    }

    */
}
 


