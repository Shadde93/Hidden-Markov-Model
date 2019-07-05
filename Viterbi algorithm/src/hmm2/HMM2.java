package hmm2;



public class HMM2 {
    
    private static int[][] specialPointsMiddle = {{1,1,1},{1,1,2},{1,2,2},{2,1,1},{2,1,2},{1,2,1},{2,2,1}, {2,2,2}};
    

    public static void main(String[] args) {
        
        for( int[] specialPointMid : specialPointsMiddle ){
            for (int i = 0; i<3;i++){
        System.out.println(specialPointMid[i]);}
}
}
}
    

