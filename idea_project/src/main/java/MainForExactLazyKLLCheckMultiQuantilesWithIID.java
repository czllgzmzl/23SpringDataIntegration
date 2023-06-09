import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Random;

/** Compare sketchByPage & sketchByChunkDivide.
 * Different ChunkSize.
*/
public class MainForExactLazyKLLCheckMultiQuantilesWithIID {
    int dataType=-233;
    static int startType=0,endType=0;
    static int pageN = 8192;
    static int N = pageN * 6712/*1000000*/, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 500; // CHECK IT
    public static int MULTI_QUANTILES = 10; // CHECK IT
    static double[] a;
//    static double[] pageMinV,pageMaxV;
    static int compaction_level;
    ArrayList<String> time_result = new ArrayList<>();
    int RESULT_LINE = 0;
    Random random = new Random(233);

    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        if (dataType == 0) {
            for (int i = 0; i < N; i++)
                a[i] = //longToResult(i&4095);
                    //Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 2) - 1) * 300);
                    random.nextGaussian();
//                    i-0.5*Math.floor(i/pageN)*pageN;
//                    i;
            return;
        }
        BufferedReader reader = null;
        if (dataType == 1||dataType==5)
            reader = new BufferedReader(new FileReader(new File("1_bitcoin.csv")));
        if (dataType == 2||dataType==6)
            reader = new BufferedReader(new FileReader(new File("2_SpacecraftThruster.txt")));
        if (dataType == 3||dataType==7)
            reader = new BufferedReader(new FileReader(new File("3_taxipredition8M.txt")));
        if (dataType == 4||dataType==8)
            reader = new BufferedReader(new FileReader(new File("4_wh.csv")));
        assert reader != null;
        reader.readLine(); // ignore first line.
        String line;
        int cntN = 0;
        while ((line = reader.readLine()) != null) {
            a[cntN++] = Double.parseDouble(line);
            if(dataType>=5)a[cntN-1]*=(1.0+random.nextGaussian()*1e-6);
            if (cntN == N) break;
        }
    }

//    public void prepareWorker(int maxSeriesByte) {
//        pageMinV = new double[pageNum];pageMaxV = new double[pageNum];
//
//        int enoughMemByte = pageN * 9;
//        for (int i = 0; i < pageNum; i++) {
//            pageMinV[i] = Long.MAX_VALUE;
//            pageMaxV[i] = Long.MIN_VALUE;
//            for (int j = 0; j < pageN; j++) {
//                double v =a[i*pageN+j];
//                pageMinV[i] = Math.min(pageMinV[i],v);
//                pageMaxV[i]=Math.max(pageMaxV[i],v);
//            }
////            worker.show();
//        }
//    }

    public int getValueActualRank(double[]sortedA,int queryN, double v){ // number of elements <= v
        int L=0,R=queryN-1;
        while(L<R){
            int mid=(L+R+1)>>>1;
            if(v<sortedA[mid])R=mid-1;
            else L=mid;
        }
        return L;
    }
    public int getValueLessThan(double[]sortedA,int queryN, double v){ // number of elements <= v
        int L=0,R=queryN-1;
        while(L<R){
            int mid=(L+R+1)>>>1;
            if(sortedA[mid]<v)L=mid;
            else R=mid-1;
        }
        return sortedA[L]<v?L:L-1;
    }
    public int getDeltaRank(double[]sortedA,int queryN, double v,int targetRank){
        int rank_L = getValueLessThan(sortedA,queryN,v)+1;
        int rank_R = getValueActualRank(sortedA,queryN,v);
//        System.out.println("\t\t\t"+targetRank+"\t\tresultLR:"+rank_L+"..."+rank_R+"\t\tresV:"+v);
        if(targetRank>=rank_L&&targetRank<=rank_R)return 0;
        else return targetRank<rank_L?(targetRank-rank_L):(targetRank-rank_R);
    }
    private long dataToLong(double data) {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }
    private double longToResult(long result) {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
    }


    public void testError(int queryN, int maxMemoryByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long full_time = 0, merge_page_time = 0;
        double avg_iteration=0;
        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }
        double ALLPrSum=0,ALLPrCount=0,FailCount=0,IterCount=0;

//        System.out.println("\t!! DataSketchesKLL_K="+DataSketchesKLL_K);
        IntArrayList rankErrList =new IntArrayList();
        for (int T = 0; T < TEST_CASE; T++) {
//            System.out.println("\tTEST_ID:"+T);
//            int L = 0, R = queryN;
            int L = LL[T], R = RR[T];
            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);
            DoubleArrayList answer = new DoubleArrayList();
            for(double deltaQ=1.0/(MULTI_QUANTILES+1),queryQ=deltaQ;queryQ<1.0-1e-9;queryQ+=deltaQ)
                answer.add(query_a[(int) (queryQ * queryN)]);

            KLLSketchLazyExact firstWorker = new KLLSketchLazyExact(maxMemoryByte);
            for (int i = L; i < R; i++) firstWorker.update(dataToLong(a[i]));
//            firstWorker.showCompact();

            int qid=0,answerApproxRankErr=0;
            for(double deltaQ=1.0/(MULTI_QUANTILES+1),queryQ=deltaQ;queryQ<1.0-1e-9;queryQ+=deltaQ) {
                int query_rank = (int) (queryQ * queryN);
                int tmpApproxRank = firstWorker.getApproxRank(dataToLong(answer.getDouble(qid)));
                int tmpApproxRank0 = firstWorker.getApproxRank(dataToLong(answer.getDouble(qid)*(1-1e-9)));
                int tmpApproxRankErr = Math.min(Math.abs(tmpApproxRank0 - query_rank),Math.abs(tmpApproxRank - query_rank));
                if(tmpApproxRank0<=query_rank&&query_rank<=tmpApproxRank)tmpApproxRankErr=0;
                qid++;
                answerApproxRankErr=Math.max(answerApproxRankErr,tmpApproxRankErr);
            }
//            System.out.println("\t\trankErr:" + answerApproxRankErr + "(" + 1.0 * answerApproxRankErr / queryN + ")");
            rankErrList.add(answerApproxRankErr);
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("CheckErrPR!\tTEST_CASE="+TEST_CASE+"\tDATASET:"+dataType+"\tqueryN:\t" + queryN+"\tmemory:\t" + maxMemoryByte+"\t\tMULTI_QUANTILES:\t"+MULTI_QUANTILES);
        KLLSketchLazyExact firstWorker = new KLLSketchLazyExact(maxMemoryByte);
        for (int i = 0; i < queryN; i++) firstWorker.update(dataToLong(a[i]));
        rankErrList.sort(Integer::compare);
        DoubleArrayList prList = new DoubleArrayList();
        for(double pr=0.5;pr<=0.90-1e-8;pr+=0.01)prList.add(pr);
        for(double pr=0.90;pr<=0.99-1e-8;pr+=0.01)prList.add(pr);
        for(double pr=0.99;pr<=0.999-1e-8;pr+=0.001)prList.add(pr);
        for(double pr=0.999;pr<=0.9999-1e-8;pr+=1e-4)prList.add(pr);
        for(double pr=0.9999;pr<=0.99999+1e-8;pr+=1e-5)prList.add(pr);
        prList.add(1.0);
//        for(double pr=0.01;pr<=0.99+1e-5;pr+=0.01)prList.add(pr);
//        for(double pr=0.991;pr<=0.999;pr+=0.001)prList.add(pr);
        firstWorker.showCompact();
        for(double pr:prList){
//            System.out.println("pr:\t"+pr+"\t\tactualPrErr:\t"+rankErrList.getInt((int)Math.floor(TEST_CASE*(pr)-1e-6))+"\t\testimated_PrErr:\t"+firstWorker.queryRankErrBound(dataToLong(0.5),pr));
            System.out.println("pr:\t"+pr+"\t\tactualMaxPrErr:\t"+rankErrList.getInt((int)Math.floor(TEST_CASE*(pr)-1e-6))+"\t\testimatedMaxPrErr:\t"+firstWorker.queryRankErrBound(Math.pow(pr,1.0/MULTI_QUANTILES))+"\tsinglePr:\t"+Math.pow(pr,1.0/MULTI_QUANTILES));
        }
        System.out.println("\n\n");

//
//        // priori
//        KLLSketchLazyEmptyForSimuCompact simulator = new KLLSketchLazyEmptyForSimuCompact(maxMemoryByte/8);
//        int[] compactNum = simulator.simulateCompactNumGivenN(queryN);
//        double tot_sig2 = simulator.getSig2()*MULTI_QUANTILES,maxErr=simulator.getMaxError();
//
//        rankErrList.sort(Integer::compare);
//        DoubleArrayList prList = new DoubleArrayList();
//        for(double pr=0.5;pr<=0.90-1e-8;pr+=0.01)prList.add(pr);
//        for(double pr=0.90;pr<=0.99-1e-8;pr+=0.01)prList.add(pr);
//        for(double pr=0.99;pr<=0.999-1e-8;pr+=0.001)prList.add(pr);
//        for(double pr=0.999;pr<=0.9999-1e-8;pr+=1e-4)prList.add(pr);
//        for(double pr=0.9999;pr<=0.99999+1e-8;pr+=1e-5)prList.add(pr);
//        prList.add(1.0);
//        for(double pr:prList){
////            System.out.println("pr:\t"+pr+"\t\tactualPrErr:\t"+rankErrList.getInt((int)Math.floor(TEST_CASE*(pr)-1e-6))+"\t\testimated_PrErr:\t"+firstWorker.queryRankErrBound(dataToLong(0.5),pr));
//            System.out.println("pr:\t"+pr+"\t\tactualMaxPrErr:\t"+rankErrList.getInt((int)Math.floor(TEST_CASE*(pr)-1e-6))+"\t\testimatedMaxPrErr:\t"+KLLSketchLazyExactPriori.queryRankErrBoundGivenParameter(tot_sig2,(int)maxErr,pr));
//        }
    }

    public void testPrioriError(int queryN, int maxMemoryByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        long full_time = 0, merge_page_time = 0;
        double avg_iteration=0;
        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }
        double ALLPrSum=0,ALLPrCount=0,FailCount=0,IterCount=0;

//        System.out.println("\t!! DataSketchesKLL_K="+DataSketchesKLL_K);
        IntArrayList rankErrList =new IntArrayList();
        for (int T = 0; T < TEST_CASE; T++) {
//            System.out.println("\tTEST_ID:"+T);
//            int L = 0, R = queryN;
            int L = LL[T], R = RR[T];
            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);
            DoubleArrayList answer = new DoubleArrayList();
            for(double deltaQ=1.0/(MULTI_QUANTILES+1),queryQ=deltaQ;queryQ<1.0-1e-9;queryQ+=deltaQ)
                answer.add(query_a[(int) (queryQ * queryN)]);

            KLLSketchLazyExactPriori firstWorker = new KLLSketchLazyExactPriori(maxMemoryByte);
            for (int i = L; i < R; i++) firstWorker.update(dataToLong(a[i]));
//            firstWorker.showCompact();

            int qid=0,answerApproxRankErr=0;
            for(double deltaQ=1.0/(MULTI_QUANTILES+1),queryQ=deltaQ;queryQ<1.0-1e-9;queryQ+=deltaQ) {
                int query_rank = (int) (queryQ * queryN);
                int tmpApproxRank = firstWorker.getApproxRank(dataToLong(answer.getDouble(qid)));
                int tmpApproxRank0 = firstWorker.getApproxRank(dataToLong(answer.getDouble(qid)*(1-1e-9)));
                int tmpApproxRankErr = Math.min(Math.abs(tmpApproxRank0 - query_rank),Math.abs(tmpApproxRank - query_rank));
                if(tmpApproxRank0<=query_rank&&query_rank<=tmpApproxRank)tmpApproxRankErr=0;
                qid++;
                answerApproxRankErr=Math.max(answerApproxRankErr,tmpApproxRankErr);
            }
//            System.out.println("\t\trankErr:" + answerApproxRankErr + "(" + 1.0 * answerApproxRankErr / queryN + ")");
            rankErrList.add(answerApproxRankErr);
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("CheckErrPR!\tTEST_CASE="+TEST_CASE+"\tDATASET:"+dataType+"\tqueryN:\t" + queryN+"\tmemory:\t" + maxMemoryByte+"\t\tMULTI_QUANTILES:\t"+MULTI_QUANTILES);

        // priori
        KLLSketchLazyEmptyForSimuCompact simulator = new KLLSketchLazyEmptyForSimuCompact(maxMemoryByte/8);
        int[] compactNum = simulator.simulateCompactNumGivenN(queryN);
        double tot_sig2 = simulator.getSig2(),maxErr=simulator.getMaxError();

        rankErrList.sort(Integer::compare);
        DoubleArrayList prList = new DoubleArrayList();
        for(double pr=0.5;pr<=0.90-1e-8;pr+=0.01)prList.add(pr);
        for(double pr=0.90;pr<=0.99-1e-8;pr+=0.01)prList.add(pr);
        for(double pr=0.99;pr<=0.999-1e-8;pr+=0.001)prList.add(pr);
        for(double pr=0.999;pr<=0.9999-1e-8;pr+=1e-4)prList.add(pr);
        for(double pr=0.9999;pr<=0.99999+1e-8;pr+=1e-5)prList.add(pr);
        prList.add(1.0);
        for(double pr:prList){
//            System.out.println("pr:\t"+pr+"\t\tactualPrErr:\t"+rankErrList.getInt((int)Math.floor(TEST_CASE*(pr)-1e-6))+"\t\testimated_PrErr:\t"+firstWorker.queryRankErrBound(dataToLong(0.5),pr));
            System.out.println("pr:\t"+pr+"\t\tactualMaxPrErr:\t"+rankErrList.getInt((int)Math.floor(TEST_CASE*(pr)-1e-6))+"\t\testimatedMaxPrErr:\t"+KLLSketchLazyExactPriori.queryRankErrBoundGivenParameter(tot_sig2,(int)maxErr,Math.pow(pr,1.0/MULTI_QUANTILES)));
        }
    }


    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForExactLazyKLLCheckMultiQuantilesWithIID main;

        main = new MainForExactLazyKLLCheckMultiQuantilesWithIID();
        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
            main.RESULT_LINE=0;
            main.prepareA(dataType);
            for (int queryN : new int[]{100000/*,200000*/})
                for (int query_mem : new int[]{/*1024*2,*/1024 * 16})
//                    main.testError(queryN, query_mem);
            main.testPrioriError(queryN, query_mem);
//                    for(double tmpPr=0.001;tmpPr>=1e-10;tmpPr/=2) {
//                        main.tmpPR=1-tmpPr;
//                        main.testError(queryN, query_mem);
//                    }
                main.RESULT_LINE++;
            System.out.println("\n-------------------------\n");
        }
//        System.out.println("\nError rate:");
//        for (String s : main.err_result)
//            System.out.println(s);
        System.out.println("\t\t\tALL_TIME:"+(new Date().getTime()-START_T));
    }
}
