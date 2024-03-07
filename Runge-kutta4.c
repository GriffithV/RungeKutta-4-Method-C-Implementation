
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// Ay'+By + G(x) = 0 

//STRUCT FOR POLYNOMINAL G(X) TERMS 
typedef struct {
    double coeff;
    double expn;
}poly_terms;

//STRUCT FOR POLYNOMINAL G(X) FUNCTION
typedef struct {
    int term_size;
    poly_terms *polyn;

}polynfunc_G;

//STRUCT FOR ODE THAT USER WANT TO CALCULATE 
typedef struct {
    
    double coefA;
    double coefB;
    polynfunc_G G_x;

}ODE;



//INITIALIZING ODE
ODE *initialize_ODE(){

    int i,n;
    //DYNAMIC MEMORY ALLOCATION FOR OUR ODE AND INITIALIZING
    ODE *d1=(ODE*)malloc(sizeof(ODE)); // 
    printf("Enter the information for ODE that you want to calculate by Runga Kutta4 \n");
    printf("The ODE has to be like ->  Ay'+ By + G(x) = 0   (G(x) is  a polynominal function)\n");
    printf("---------------------------------------------------------------------------------------\n");
    printf("Please enter coefficent for A :");   scanf("%lf",&d1->coefA);
    printf("\nPlease enter coefficent for B :");   scanf("%lf",&d1->coefB);
    
    printf("\nPlease enter the number of terms for polynomial func G(X) :");    scanf("%d",&n);
    d1->G_x.term_size=n;
    d1->G_x.polyn=malloc(n*sizeof(poly_terms)); // DYNAMIC MEMORY ALLOCATION FOR G(X) TERMS
    printf("\n");
     //INITIALIZING EACH TERMS OF G(X)
    for(i=0;i<n;i++){ 
        printf("\nCoefficent for %d.Term  Ax^B enter A:",i+1); scanf("%lf",&((d1->G_x.polyn+i)->coeff));
        printf("\nExponential for %d.Term  Ax^B   enter B:",i+1); scanf("%lf",&((d1->G_x.polyn+i)->expn));
    }

    return d1;
}



//PRINTING THE GIVEN ODE AS OUTPUT
void print_given_ODE(ODE *d){
    int i;
    printf("\nThe Given ODE is : %.2lfy' + %.2lfy ",d->coefA,d->coefB);
    for(i=0;i<d->G_x.term_size;i++){
        printf("+ %.2lfx^%.2lf",(d->G_x.polyn+i)->coeff,(d->G_x.polyn+i)->expn);
 
    }
    printf(" = 0\n");
}


//CALCULATING Y(x) FOR X VALUE WHICH USER WANTS BY RUNGE KUTTA 4 
void runge_kutta_4(ODE *d){   
    double actual_value,loss,x0,y0,result_x,h,k1,k2,k3,k4;

    printf("\nEnter the initial values for f(x)=y\n x = "); scanf("%lf",&x0); printf(" y = "); scanf("%lf",&y0);
    printf (" y(%.0lf)=%.0lf\n",x0,y0);
    printf("\nTo calculate f(x) , enter a 'x' value : "); scanf("%lf",&result_x);
    printf ("Enter the step size : "); scanf("%lf",&h);
    printf ("To estimating errors , enter the actual value for y(%.2lf) = ",result_x); scanf("%lf",&actual_value);
    printf ("Calculating y(%.1lf) by Range Kutta4 with  step size -> h = %.1lf\n",result_x,h);
    
    int i,j;
    double polysum=0;
    double iter= ((result_x-x0)/h); // number of steps
    for(i=0;i<iter;i++){
        for(j=0;j<d->G_x.term_size;j++){
            polysum+=((d->G_x.polyn+j)->coeff)*pow(x0,(d->G_x.polyn+j)->expn); //Calculates G(xi) value by using struct of polynominal G(X)
        }
        //CALCULATING K1
        k1=((y0*-1.0*(d->coefB)) + (-1.0*polysum))/d->coefA;
        polysum=0;
        for(j=0;j<d->G_x.term_size;j++){
            polysum+=((d->G_x.polyn+j)->coeff)*pow((x0+0.5*h),(d->G_x.polyn+j)->expn);
        }
        //CALCULATING K2
        k2=(((y0+0.5*k1*h)*-1.0*(d->coefB)) + (-1.0*polysum))/d->coefA;
        polysum=0;
        for(j=0;j<d->G_x.term_size;j++){
            polysum+=((d->G_x.polyn+j)->coeff)*pow((x0+0.5*h),(d->G_x.polyn+j)->expn);
        }
        //CALCULATING K3
        k3=(((y0+0.5*k2*h)*-1.0*(d->coefB)) + (-1.0*polysum))/d->coefA;
        polysum=0;
        for(j=0;j<d->G_x.term_size;j++){
            polysum+=((d->G_x.polyn+j)->coeff)*pow((x0+h),(d->G_x.polyn+j)->expn);
        }
        //CALCULATING K4
        k4=(((y0+k3*h)*-1.0*(d->coefB)) + (-1.0*polysum))/d->coefA;
        polysum=0;
        //CALCUTAING  y(xi + h)
        y0=y0+((k1+(2*k2)+(2*k3)+k4)*h)/6; 
        x0+=h;
        loss=actual_value-y0;//CALCULATING ERRORS FOR EACH ITERATION
        if(loss<0) loss=loss*-1; //Absolute Error 

        printf ("\n %d iter -> y(%.2lf)=%.13lf  estimated error : %.9lf\n",i+1,x0,y0,loss);
    }

    printf("\n Approximately , the result for y(%.2lf)=  %lf",result_x,y0); // prints final result 
    printf("\n with error :  %.13lf ",loss); // with error 


}


int main(){
   //INITIALIZING ODE AND CALCULATING Y(X) by RUNGE KUTTA4
   ODE *d=initialize_ODE();
   print_given_ODE(d);
   runge_kutta_4(d);

   //DEALLOCATING MEMORY
   free(d->G_x.polyn);
   free(d);

  return 0;
}
