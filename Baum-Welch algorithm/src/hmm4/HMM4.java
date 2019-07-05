package hmm4;

import java.util.*;
import java.io.*;
import java.math.*;

public class HMM4 {
                   
    private int iteration = 0;
    private double [][] a;
    private List<Integer> dimA;
    private int numStates;
    private double [][] b;
    private List<Integer> dimB;
    private double [] pi;
    private List<Integer> dimPI;
    private List<Integer> finalDim;
    private int[] observations;
    private int numObs;
    private double[][][] digamma;
    private double[][] gamma;
    private double[]scaleList;
    public List<Integer> getdimA(){
        return this.dimA;
    }
    public List<Integer> getdimB(){
	return this.dimB;
    }
    public List<Integer> getdimPI(){
        return this.dimPI;
    }
    public List<Integer> getfinalDim(){
        return this.finalDim;
        }
    public double[][] getA(){
	return this.a;
	}
    public double[][] getB(){
	return this.b;
	}
    public double[] getPI(){
	return this.pi;
	}
    public int[] getObs(){
	return this.observations;
    }
    public int getnumStates(){
	return this.numStates;
    }
    public double[][] getGamma(){
	return this.gamma;
    }
    public double[][][] getDiGamma(){
	return this.digamma;
    }
 
    public void HMM4(){
        this.a=a;
        this.b=b;
        this.dimA=dimA;
        this.dimB=dimB;
        this.dimPI=dimPI;
        this.finalDim = finalDim;
        this.numObs = numObs;
        this.numStates = numStates;
        this.observations = observations;
        this.pi = pi;        
        this. digamma = digamma;
        this.gamma = gamma;
        this. scaleList = scaleList;
                    
    }
                
    private void setA(List<Double> val, List<Integer> dim) {
        this.numStates = dim.get(0);
        this.a = new double[numStates][numStates];
        this.dimA = new ArrayList<> (dim);
        int k = 0;

        for(int i = 0; i < numStates; i++) {
            for(int j = 0; j < numStates; j++) {
                // read information from somewhere
                this.a[i][j] = val.get(k);
                k++;
            }
        }
    		//System.out.println(" A has been set");
    }
    

    private void setB(List<Double> val, List<Integer> dim) {
        this.b = new  double[dim.get(0)][dim.get(1)];
        this.dimB = new ArrayList<Integer> (dim);
        int k = 0;

	for(int i = 0; i < dim.get(0); i++) {
    		for(int j = 0; j < dim.get(1); j++) {
                // read information from somewhere
                this.b[i][j] = val.get(k);
    		k++;
    		}
    	}
    		//System.out.println("B has been set");
    }

    private void setPI(List<Double> val, List<Integer> dim) {
   	this.dimPI = new ArrayList<Integer> (dim);
	this.pi = new double[dim.get(1)];
            
            for(int i = 0; i < dim.get(1); i++) {
    		this.pi[i] = val.get(i);
            }
            //System.out.println(" PI has been set");
    }

    private void setObs(int l, List<Integer> val){
    	this.observations = new int[l];
    	this.numObs = l;

	for(int i = 0; i < l; i++) {
            this.observations[i] = val.get(i);
    	}
            //System.out.println("Observation sequence set");
    }
    
    public void readFile(){
	List<Double> array = new ArrayList<Double>();
        List<Integer> dimensions = new ArrayList<Integer>();
	List<Integer> observations = new ArrayList<Integer>();
	int numObs = 0;
	int numLine =0;
	int n = 0;
	try {
            Scanner numFile = new Scanner(new File("hmm4_01.in"));
            //Scanner numFile = new Scanner(System.in);
            while (numFile.hasNextLine()){
                String line = numFile.nextLine();
    		Scanner sc = new Scanner(line);
                sc.useDelimiter(" ");
		if (numLine<3){
    			//System.out.println(line);
                    while(sc.hasNext()) {
    			if (n < 2 ) {
                            int i = Integer.parseInt(sc.next());
                            dimensions.add(i);
                            n++;
                        }
       			else {
                            Double d = Double.parseDouble(sc.next());
                            //System.out.println(d);
                            array.add(d);
                            n++;
        		}
                    }
                    switch (numLine) {            
                        case 0 : setA(array, dimensions); 
                            array.clear();
                            dimensions.clear();
                            n=0;
                                break;
                        case 1 : setB(array, dimensions);
                            array.clear();
                            dimensions.clear();
                            n=0;
                		 break;
                        case 2 : setPI(array, dimensions);
                            array.clear();
                            dimensions.clear();
                            n=0;  
                		 break;
                    }
                    sc.close();
                    numLine++;
    		}
                else {
                    while(sc.hasNext()) {
                        if (n < 1 ) {
                            numObs = Integer.parseInt(sc.next());
                            n++;
                        }
       			else {
                            observations.add(Integer.parseInt(sc.next()));
        		}
                    }
                    n=0;
                    numLine=0;
                    setObs(numObs, observations);
                }
                //System.out.print(array);
            }
        }
	catch (Exception e) {
        }
    }

    public double[] calculate(List<Integer> dimA, List<Integer> dimB, double [] a, double [][] b){
        double c = 0;
        double rowsum = 0;
        double prod[]  = new double[dimB.get(1)];
        for (int k = 0; k < dimB.get(1); k++){
		for (int r =0; r < dimB.get(0); r++){
                    c = a[r] * b[r][k];
                    //System.out.println(c);
                    prod[k] = prod[k] + c;
		}
		c = 0;
	}
	List<Integer> newDim = new ArrayList<Integer>();
	newDim.add(dimA.get(0));
	newDim.add(dimB.get(1));
	this.finalDim = new ArrayList<Integer> (newDim);
            return prod;
    }         
                
    public void BaumWelch(double[][] a, double[][] b, double[] pi, int[] o){

    	this.a = a;
    	this.b = b;
    	this.pi = pi;
    	double[][] alphaMatrix = alphaPass(o);
		double[][] betaMatrix = betaPass(o);
		gamma(alphaMatrix, betaMatrix, o);
		double[][] newA = estimateA(gamma, digamma);
		double[][] newB = estimateB(gamma, digamma);

		double newLogProb = 0;
		for (int t=0; t< o.length; t++) {
			newLogProb += Math.log(scaleList[t]);
		}
		newLogProb = 0 - newLogProb;

		//Check if iterate or completed
		iteration++;
		int limit = 500;
		double logProb = Double.NEGATIVE_INFINITY;

		if (iteration >= limit || newLogProb < logProb) {
                    printAnswer(newA, newB);
		} else {
                    BaumWelch(newA, newB, pi, o);
		}
    }

    	private void printAnswer(double[][] a , double[][] b){
            String finalanswerString = "";
            finalanswerString = (dimA.get(0) + " " + dimA.get(0) + " ");
			//print final estimate
            for(int j=0; j< numStates; j++){
                for (int i = 0; i < numStates; i++){
	  	finalanswerString += (a[j][i] + " ");
	    }
        }

	      finalanswerString += ("\n" + dimB.get(0) + " " +  dimB.get(1) + " ");


	     for(int j=0; j< numStates; j++){

	     for (int i = 0; i < dimB.get(1); i++){
	  	  finalanswerString += (b[j][i] + " ");
	      }
	      }


	  System.out.println(finalanswerString);

    	}


            //taken from stamp
  public double[][] alphaPass(int[] o) {
            
            //
	    int T = o.length;
	    double[][] alphas = new double[numStates][T];
	    double[] scaleList = new double[T];
	    scaleList[0] = 0;
	    /* initialization (time 0) */
	    for (int i = 0; i < numStates; i++){
	      alphas[i][0] = pi[i] * b[i][o[0]];
	      scaleList[0] += alphas[i][0];

	      }

	    //scaling
	    scaleList[0] = 1/scaleList[0];
	    for (int i =0; i<numStates; i++) {
	    	alphas[i][0] = alphas[i][0]*scaleList[0];
            }

	    /* induction */
	    for (int t = 1; t < T; t++) {
	    	scaleList[t] = 0;
	      for (int i = 0; i < numStates; i++) {
				alphas[i][t] = 0;
				
				for (int j = 0; j < numStates; j++){
		  		alphas[i][t] += (alphas[j][t-1] * a[j][i]);
		  		}
		  	alphas[i][t] *= b[i][o[t]];
		  	scaleList[t] +=	alphas[i][t];
	      }
	     //scaling

	    scaleList[t] = 1/scaleList[t];

	    for (int i=0; i<numStates; i++) {
	    	alphas[i][t] = scaleList[t]*alphas[i][t];
	    }

	    }

	    //Prints all alphas at time t
	     //System.out.println("Alpha starts here");
	     /*for(int j=0; j<8; j++){
	     System.out.print("\n");
	     for (int i = 0; i < numStates; i++){
	  	  System.out.println(alphas[i][j]);
	      }
	      }*/
	      this.scaleList = scaleList;
	      return alphas;
}
  //taken from stamp
  public double[][] betaPass(int[] o) {
	    int T = o.length;
	    double[][] betas = new double[numStates][T];
	        
	    /* initialization (time 0) */
	    for (int i = 0; i < numStates; i++){
	      betas[i][T-1] = scaleList[T-1];
	    }


	    /* induction */
	    for (int t = T - 2; t >= 0; t--) {
	      for (int i = 0; i < numStates; i++) {
			betas[i][t] = 0;
			for (int j = 0; j < numStates; j++){
		  	betas[i][t] += (betas[j][t+1] * a[i][j] * b[j][o[t+1]]);
		  	}
		  	betas[i][t] = betas[i][t] * scaleList[t];
	      	}
	      	
	    	}


	   // Prints all betas at time t
	    /*System.out.println("Beta starts here");
	     for(int j=0; j< o.length; j++){
	     System.out.print("\n");
	     for (int i = 0; i < numStates; i++){
	  	  System.out.println(betas[i][j]);
	      }
	      }*/


	    	return betas;
	  	}
//taken from stamp
  public void gamma(double[][] alpha, double[][] beta, int [] o){
  	int T = o.length;
  	double[][][] digamma = new double[numStates][numStates][T-1];
  	double[][] gamma = new double[numStates][T];
  	double denom = 0;


        //
  	for (int t=0; t < T-1; t++){
  		denom = 0;
  		for (int i = 0; i < numStates; i++) {
  			for (int j = 0; j < numStates; j++) {
  				denom += (alpha[i][t] * a[i][j] * b[j][o[t+1]] * beta[j][t+1]);
  			}
  
  			}
  			
  		for (int i = 0; i < numStates ;i++ ) {
  			gamma[i][t] = 0;
  			for (int j=0; j< numStates;j++ ) {
  				digamma[i][j][t] = (alpha[i][t]*a[i][j] * b[j][o[t+1]] * beta[j][t+1])/denom;
  				gamma[i][t] += digamma[i][j][t];
  				
  			}
	
  		}
  		}

  		denom = 0;
  		for (int i = 0; i < numStates; i++) {
  			denom += alpha[i][T-1];
  			
  		}
  		for (int i = 0; i < numStates; i++) {
  	
  			gamma[i][T-1] = alpha[i][T-1]/denom;
  		}

  		this.gamma = gamma;
  		this.digamma = digamma;
  	
  		   // Prints all digamma at time t
  	/*System.out.print("gamma starts here");
	     for(int j=0; j< T; j++){
	     System.out.print("\n");
	     for (int i = 0; i < numStates; i++){
	  	  System.out.println(gamma[i][j]);
	      }
	      }
	  System.out.println();*/
}
  
  //taken from stamp
 public double[][] estimateA(double[][]gamma, double[][][]digamma){
	double[][] a = new double[numStates][numStates];
	double val = 0;
	double div = 0;

		for(int i = 0; i < numStates; i++) {
    		for(int j = 0; j < numStates; j++) {
    			val = 0;
    			div = 0;//
    			for (int t=0; t < numObs-1;t++) {
    				val += digamma[i][j][t];
    				div += gamma[i][t];
    				
    			}

    			a[i][j] = val/div;
    			
    			}
    		}


    /*System.out.print("A starts here");
	     for(int j=0; j< numStates; j++){
	     System.out.print("\n");
	     for (int i = 0; i < numStates; i++){
	  	  System.out.println(a[i][j]);
	      }
	      }*/



	return a;
}
 //taken from stamp
 public double[][] estimateB(double[][]gamma, double[][][]digamma){
	double[][] b = new double[numStates][dimB.get(1)];
	double val = 0;
	double div = 0;
            //
		for(int i = 0; i < numStates; i++) {
    		for(int j = 0; j < dimB.get(1); j++) {
    			val = 0;
    			div = 0;
    			for (int t=0; t < numObs;t++) {
    				if(observations[t] == j){
    				 	val += gamma[i][t];
    				}
    				div += gamma[i][t];

    			
    				
    			}

    			b[i][j] = val/div;
    			
    			}
    		}


    /*System.out.print("B starts here");
	     for(int j=0; j< dimB.get(1); j++){
	     System.out.print("\n");
	     for (int i = 0; i < numStates; i++){
	  	  System.out.println(b[i][j]);
	      }
	      }*/



	return b;
}
 //taken from stamp
public double [] estimatePI (double[][] gamma){
	double[] pi = new double[numStates];
	for (int i = 0; i < numStates ;i++ ) {
		pi[i] = gamma[i][0];
		
	}
	this.pi = pi;
        
        	     /*for(int j=0; j< numStates; j++){
	     System.out.print("\n");
	  	  System.out.println(pi[j]);
	      }*/
	return pi;
}

  



		public static void main(String args[]) {

				

				HMM4 h = new HMM4();
				h.readFile();
				double[][] alphaMatrix = h.alphaPass(h.getObs());
				double[][] betaMatrix = h.betaPass(h.getObs());
                                
				//double[][] digamma = h.digamma(alphaMatrix, betaMatrix, h.getObs());
				h.gamma(alphaMatrix, betaMatrix, h.getObs());
                                double[] pi = h.estimatePI(h.gamma);

				double[][] newA = h.estimateA(h.getGamma(), h.getDiGamma());
				double[][] newB = h.estimateB(h.getGamma(), h.getDiGamma());
	      		//double [][] y = h.betaPass(h.getObs());
				//System.out.println(x);
                                

				
				h.BaumWelch(h.getA(), h.getB(), h.getPI(), h.getObs());


		}

}
