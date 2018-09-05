#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <tuple>
#include <algorithm>


using namespace std;


vector<double> linspace(double start, double end, double steps){
    vector<double> x;
    double step_length = (end-start)/(steps);
    for (double i = start; i < end+1; i+=step_length){
        x.push_back(i);
    }
    return x;
}

double f(double x){
    return 100*exp(-10*x);
}

double maks(vector<double> y){
    double storst;
    for (int i = 1; i<y.size()-1; i+=1){
        if (y[i] > storst){
            storst = y[i];
        }
    }
    return storst;
}

//b)
tuple<vector<double>,vector<double> ,vector<double>> task_b(int nlim){
    vector<double> a(nlim-1, 0.0);
    vector<double> b(nlim, 0.0);
    vector<double> c(nlim-1, 0.0);    

    b[0] = rand() % 100;
    for(int i=0;i<(nlim-1);i+=1){
        a[i] = rand() % 100;
        b[i+1] = rand() % 100;
        c[i] = rand() % 100;
    }
    return {a, b, c};
}


//c)
tuple<vector<double>,vector<double> ,vector<double>> task_c(int nlim){
    vector<double> a(nlim-1, 0.0);
    vector<double> b(nlim, 0.0);
    vector<double> c(nlim-1, 0.0);    

    b[0] = 2;
    for(int i=0;i<(nlim-1);i+=1){
        a[i] = -1;
        b[i+1] =  2;
        c[i] = -1;
    }
    return {a, b, c};
}

int main(int argc,char* argv[]){
    ofstream myfile;
    
    float n = pow(10.0, atof(argv[1])); 
    float h = 1.0/(n+1);
    
    double s;
    

    vector<double> x = linspace(0,1,n-1);
    vector<double> v(n, 0.0);
    vector<double> b_tilde(n, 0.0);
    vector<double> u(n, 0.0);
    vector<double> eps(n, 0.0);

    u[0] = pow(10, -15);
    
    for(double i = 0; i<n; i+=1){
        b_tilde[i] = pow(h, 2.0)*f(x[i]);
               
    }

//installing vectors a, b and c for task b or c
    auto tas = task_c(n);
    vector<double> a = get<0>(tas);
    vector<double> b = get<1>(tas);
    vector<double> c = get<2>(tas);

 
    
//forward substitution
    for(int i = 1; i<n; i+=1){
        s = a[i-1]/b[i-1];  
        b[i] = b[i]-s*c[i-1];
        b_tilde[i]=b_tilde[i]-s*b_tilde[i-1];
        
    } 

//backward substitution
    v[n-1]=b_tilde[n-1]/b[n-1];
    
    for(int i = 1; i<n; i+=1){
        int j = n-1-i;
        v[j]=(b_tilde[j]-c[j]*v[j+1])/b[j];
        
    }

//exact solution
    for(int i = 1; i<n; i+=1){
        u[i]=(1-(1-exp(-10))*x[i]-exp(-10*x[i]));

    }
   
    myfile.open ("task.txt");
    
    for(int i = 0; i<n; i+=1){
        myfile << x[i] << " " << v[i] << " " << u[i] << endl;
    }
    myfile.close();

    for(int i =0; i<n; i+=1){
        eps[i] = log10(abs((v[i]-u[i])/u[i]));

    }
    
    double max = maks(eps);
    myfile.open("taskdn"+string (argv[1])+".txt");
    
    myfile << max << " " << h;

    
    

    return 0;
}   
