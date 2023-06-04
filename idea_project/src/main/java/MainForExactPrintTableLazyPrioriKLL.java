import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Date;


public class MainForExactPrintTableLazyPrioriKLL {
    public static int Pr_NUM = 64;
    static double errST = 0.5, errED = 5e-4, nST = 1e5, nED = 1e9, nNum = 65, mST = 8 * 1024, mED = 1024 * 1024, mNum = 64;
    static DoubleArrayList prList = null;

    static DecimalFormat dfPr = new DecimalFormat("#.##E0");
    static DecimalFormat dfPass = new DecimalFormat("0.00");
    static DecimalFormat dfM = new DecimalFormat("#.##");

    public double[] testError(int queryN, int maxMemoryByte/*,int resultID*/) throws IOException {
//        System.out.println("FindFixPR!\tTEST_CASE=" + TEST_CASE + "\tDATASET:" + dataType + "\tqueryN:\t" + queryN + "\tmemory:\t" + maxMemoryByte);
//        System.out.println("\t\tPr\testiF(N,Mem,Pr)");
////
//        double avgF = 0,avgPass=0;
//        double[] query_a = new double[queryN];
//        simuWorker=new KLLSketchLazyEmptyForSimuCompact(maxMemoryByte/8);
//        for(double fixPr:prList) {
//            double estiPass = PrioriBestPrHelper.evaluatePrFromScratch(maxMemoryByte/8,fixPr,queryN,merge_buffer_ratio,multi_quantile);
//            System.out.println("\t\t" + fixPr+ "\t" + estiPass );
//        }
        PrioriBestPrHelper helper = new PrioriBestPrHelper(maxMemoryByte / 8, queryN, 0, 1);
        double[] findPrResult = helper.findBestPr(1e-6, 1e-6, 5e-1, queryN);
//        if(resultID==0)// bestPr
//        System.out.print(dfPr.format(1-findPrResult[0])+"\t");
//        else System.out.print(dfPass.format(findPrResult[1])+"\t");
        return new double[]{findPrResult[0], findPrResult[1]};
    }

    public static boolean noNeed(int i,int j){
//        return (i < mNum / 4 && j >= nNum / 4 * 3 + i) || (i >= mNum / 4 * 3 && j < (i - mNum / 4 * 3 + 1));
        return !((i < mNum / 4 && j >= nNum / 4 * 3 + i) || (i >= mNum / 4 * 3 && j < (i - mNum / 4 * 3 + 1)));
    }
    public static void main(String[] args) throws IOException {
        prList = new DoubleArrayList();

        for (double i = 0, pNum = Pr_NUM, err = errST, rela = Math.exp(Math.log(errED / errST) / (pNum - 1)); i < pNum; i++, err *= rela) {
            prList.add(1 - err);
//            System.out.println("\t\t\t+\t+delta:\t"+err+"\t\t"+Math.log(err/errED)*(pNum-1)/Math.log(errST/errED)+"\t\t"+Math.log(errST/errED));
        }

        long START_T = new Date().getTime();
        MainForExactPrintTableLazyPrioriKLL main;

        main = new MainForExactPrintTableLazyPrioriKLL();


        IntArrayList NN = new IntArrayList();
        for (double i = 0, n = nST, rela = Math.exp(Math.log(nED / nST) / (nNum - 1)); i < nNum; i++, n *= rela)
            NN.add((int) Math.round(n));
        System.out.println(NN);
        IntArrayList MM = new IntArrayList();
        for (double i = 0, m = mST, rela = Math.exp(Math.log(mED / mST) / (mNum - 1)); i < mNum; i++, m *= rela)
            MM.add((int) Math.round(m));
        System.out.println(MM);

        double[][] MNPr = new double[(int) mNum][(int) nNum];
        double[][] MNPass = new double[(int) mNum][(int) nNum];

        for (int i = 0; i < mNum; i++)
            for (int j = 0; j < nNum; j++)
                if (noNeed(i,j))
                {MNPr[i][j] = 1;MNPass[i][j] = 0;}
                else {
                    double[] tmp = main.testError(NN.getInt(j), MM.getInt(i));
                    MNPr[i][j] = tmp[0];
                    MNPass[i][j] = tmp[1];
                    System.out.println("\t\t" + MM.getInt(i) + "\t" + NN.getInt(j) + "\t\t" + Arrays.toString(tmp));
                }


        System.out.println("X: M\tY: N");
        System.out.print("\t\t\tM(KB):\t\t");
        for (int queryMem : MM) System.out.print(dfM.format(queryMem / 1024.0) + "\t");
        System.out.println();
        System.out.print("\t\t\tN:\t\t");
        for (int queryN : NN) System.out.print(queryN + "\t");
        System.out.println();

        System.out.println("M-N-Pass");
        System.out.print("\tM(KB)\\\t\t\tN\t");
        for (int queryN : NN) System.out.print(queryN + "\t");
        System.out.println();
        for (int i = 0; i < mNum; i++) {
            System.out.print("\t" + dfM.format(MM.getInt(i) / 1024.0) + "\t" + "\t\t\t");
            for (int j = 0; j < nNum; j++)
                System.out.print(dfPass.format(MNPass[i][j]) + "\t");
            System.out.println();
        }
        System.out.println();
        System.out.println("M-N-δ");
        System.out.print("\tM(KB)\\\t\t\tN\t");
        for (int queryN : NN) System.out.print(queryN + "\t");
        System.out.println();
        for (int i = 0; i < mNum; i++) {
            System.out.print("\t" + dfM.format(MM.getInt(i) / 1024.0) + "\t" +  "\t\t\t");
            for (int j = 0; j < nNum; j++)
                System.out.print(dfPr.format(1.0 - MNPr[i][j]) + "\t");
            System.out.println();
        }

//            System.out.println("M-N-Pass");
//            for(int i=0;i<mNum;i++) {
//                int queryMem=MM.getInt(i);
//                System.out.print("\t"+dfM.format(queryMem/1024.0)+"\t"+i+"\t\t");
//                for(int j=0;j<nNum;j++) {
//                    if((i<mNum/4&&j>=nNum/4*3+i)||(i>=mNum/4*3&&j<(i-mNum/4*3+1)))
//                        System.out.print(dfPass.format(0)+"\t");
//                    else
//                        main.testError(NN.getInt(j), queryMem, 1);
//                }
//                System.out.println();
//            }
//            System.out.println();
//            System.out.println();
//            System.out.println("M-N-Pr");
//            for(int i=0;i<mNum;i++) {
//                int queryMem=MM.getInt(i);
//                System.out.print("\t"+dfM.format(queryMem/1024.0)+"\t"+i+"\t\t");
//                for(int j=0;j<nNum;j++) {
//                    if((i<mNum/4&&j>=nNum/4*3+i)||(i>=mNum/4*3&&j<(i-mNum/4*3+1)))
//                        System.out.print(dfPr.format(0)+"\t");
//                    else
//                        main.testError(NN.getInt(j), queryMem, 0);
//                }
//                System.out.println();
//            }
        System.out.println("\t\t\tALL_TIME:" + (new Date().getTime() - START_T));
    }
}
