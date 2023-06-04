import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Date;
import java.util.Random;


public class MainForExactPrintTableGK {
    int dataType = -233;
    static int startType = 0, endType = 0;
    double fixPr;
    static int pageN = 8192;
    public static int TEST_CASE = 1,QUANTILE_PER_TEST=1;
    static double errST = 0.5, errED = 5e-4, nST = 1e5, nED = 1e9, nNum = 16*4+1, mST = 8 * 1024, mED = 1024 * 1024, mNum = 64;
    static int N = (int)nED/3*4;
    static double[] a;
    //    static double[] pageMinV,pageMaxV;
    static int compaction_level;
    ObjectArrayList<StringBuilder>RESULT=new ObjectArrayList<>();
    int RESULT_LINE = 0;
    Random random = new Random(233);
    static DecimalFormat dfPass = new DecimalFormat("0.00");
    static DecimalFormat dfM = new DecimalFormat("#.##");


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
        if (dataType == 4) {
            for (int i = 0; i < N; i++)
                a[i] = i;
            return;
        }
        BufferedReader reader = null;
        if (dataType == 1)
            reader = new BufferedReader(new FileReader(new File("1_bitcoin.csv")));
        if (dataType == 2)
            reader = new BufferedReader(new FileReader(new File("2_SpacecraftThruster.txt")));
        if (dataType == 3)
            reader = new BufferedReader(new FileReader(new File("3_taxipredition8M.txt")));
        if (dataType == 4)
            reader = new BufferedReader(new FileReader(new File("4_wh.csv")));
        assert reader != null;
        reader.readLine(); // ignore first line.
        String line;
        int cntN = 0;
        while ((line = reader.readLine()) != null) {
            a[cntN++] = Double.parseDouble(line);
            if (cntN == N) break;
        }
    }
    private long dataToLong(double data) {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }
    private double longToResult(long result) {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
    }


    public double testError(int queryN, int maxMemoryByte) throws IOException {
        System.out.println("\t\t\ttestError\t"+queryN+"\t"+maxMemoryByte/1024+"KB");
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

        for (int T = 0; T < TEST_CASE; T++) {
            int L=LL[T],R=RR[T];

            merge_page_time -= new Date().getTime();
            GKSketchForExact lazyWorker = new GKSketchForExact(maxMemoryByte);

            for (int i = L; i < R; i++)
                lazyWorker.insert(a[i]);
//            lazyWorker.show();
//            lazyWorker.showCompact();

//            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
//            Arrays.sort(query_a);

            int q_count=QUANTILE_PER_TEST;
            double ratioPerQuery = 1.0 / (q_count * TEST_CASE);
            for(int INNER_T=0;INNER_T<q_count;INNER_T++)
            {
                double q = random.nextDouble();
                int query_rank1 = (int) Math.floor(q * (queryN-1)+1),query_rank2 = (int) Math.ceil(q * (queryN-1)+1);
                long last_n=queryN;

                double[] iterate_result;// = tDigestWorker.findResultRange(query_rank1,query_rank2);
                iterate_result =lazyWorker.getFilter(0,0,0,0,query_rank1,query_rank2);
                last_n = lazyWorker.findMaxNumberInRange(iterate_result[0], iterate_result[1]);
//                System.out.println("?????2");
                avg_iteration+=ratioPerQuery;
                int MMP=0;
                while(iterate_result[0]<iterate_result[1]&&iterate_result.length==2) {
                    if(++MMP>20)break;
                    avg_iteration+=ratioPerQuery;
                    GKSketchForExact cntWorker = new GKSketchForExact(maxMemoryByte);

                    double valL=iterate_result[0],valR=iterate_result[1];
                    int CountOfLessThanValL = 0,CountOfValL=0,CountOfValR=0;
//                    System.out.println("\t(start GK Iteration:"+"\t\tvalL,R:"+valL+","+valR+"\t\t\tMMP:\t"+MMP+"\t\tn:"+last_n);
                    if(last_n<=maxMemoryByte/8){
//                        System.out.println("\t\t\tlastPass\t\tmemUsedRate:\t"+1.0*last_n*8/maxMemoryByte+"\t\tlast_n:\t"+last_n);
                        DoubleArrayList tmpList = new DoubleArrayList((int)last_n);
                        for(int i=L;i<R;i++) {
                            if (a[i] >= valL && a[i] <= valR) tmpList.add(a[i]);
                            else if (a[i] < valL) CountOfLessThanValL++;
                        }
                        if(tmpList.size()>maxMemoryByte/8)System.out.println("!!!!!!!FAILED ALL VALUE TOO MUCH");
                        int k1=query_rank1-CountOfLessThanValL-1,k2=query_rank2-CountOfLessThanValL-1;
//                        System.out.println("\t\t\t\t use brute force. k1,k2:\t"+k1+","+k2+"\t\tCountOfLessThanValL:"+CountOfLessThanValL+"\t\tquery_rank:"+query_rank1+","+query_rank2);
                        tmpList.sort(Double::compare);
                        iterate_result[0] = tmpList.getDouble(k1);
                        iterate_result[1] = tmpList.getDouble(k2);
                        break;
                    }else {
                        for (int i = L; i < R; i++) {
                            if (a[i] > valL && a[i] < valR) cntWorker.insert(a[i]);
                            else if (a[i] < valL) CountOfLessThanValL++;
                            else if(a[i]==valL)CountOfValL++;
                            else if(a[i]==valR)CountOfValR++;
                        }
//                        cntWorker.show();
//                        System.out.println("\t\t\t\t\t\tCountOfLessThanValL:" + CountOfLessThanValL + "\t\tcntSketch_N:" + cntWorker.totN);
                        int cntRank1 = query_rank1 - CountOfLessThanValL;
                        int cntRank2 = query_rank2 - CountOfLessThanValL;
//                        System.out.println("\t\t\t\t\t\tcntRank:" + cntRank1 + " " + cntRank2);
//                        iterate_result = cntWorker.findResultRange(cntRank1, cntRank2);
                        iterate_result =cntWorker.getFilter(CountOfValL,CountOfValR,valL,valR,cntRank1,cntRank2);
                        iterate_result[0] = Math.max(iterate_result[0], valL);
                        iterate_result[1] = Math.min(iterate_result[1], valR);
                        last_n = cntWorker.findMaxNumberInRange(iterate_result[0], iterate_result[1]);
//                        System.out.println("\t\t\t\tcntL,R:" + iterate_result[0] + "," + iterate_result[1] + "\t\t\tCountOfLessThanValL:" + CountOfLessThanValL + "\t\tcntRank:" + cntRank1 + " " + cntRank2+"\t\tlastN:"+last_n);
                    }
                }
                double exact_quantile_v = (iterate_result[0]+iterate_result[1])*0.5;
//                System.out.println("\t\tq:"+q+"\texact_quantile="+exact_quantile_v);
            }
        }
        System.out.println("\t\t\t\t\t"+avg_iteration);
        return avg_iteration;
    }

    public double[] testErrorBatch(int queryN, IntArrayList MM) throws IOException {
        DoubleArrayList ans=new DoubleArrayList();
        for(int m:MM)
            ans.add(testError(queryN,m));
        return ans.toDoubleArray();
    }

//    public double[] testErrorBatch(int queryN, IntArrayList MM) throws IOException {
//        System.out.println("\t\tcntN:\t"+queryN);
//        int BATCH = MM.size();
//        long full_time = 0, merge_page_time = 0;
//        double[] avg_iteration=new double[BATCH];
////        for(int i=0;i<BATCH;i++)avg_iteration[i]=-1;if(queryN>=1)return avg_iteration;
//
//        int[] LL = new int[TEST_CASE];
//        int[] RR = new int[TEST_CASE];
//        Random random = new Random(233);
//        for (int i = 0; i < TEST_CASE; i++) {
//            LL[i] = random.nextInt(N - queryN + 1);
//            RR[i] = LL[i] + queryN;
//        }
//
//        for (int T = 0; T < TEST_CASE; T++) {
//            int L=LL[T],R=RR[T];
//
//            merge_page_time -= new Date().getTime();
//            ObjectArrayList<MRLSketchLazy> lazyWorkers = new ObjectArrayList<>(BATCH);
//            for(int m:MM)lazyWorkers.add(new MRLSketchLazy(m));
////            MRLSketchLazy lazyWorker = new MRLSketchLazy(maxMemoryByte);
//
//            for (int i = L; i < R; i++)
//                for(MRLSketchLazy lazyWorker:lazyWorkers)
//                    lazyWorker.update(dataToLong(a[i]));
//
//            int q_count=QUANTILE_PER_TEST;
//            double ratioPerQuery = 1.0 / (q_count * TEST_CASE);
//            for(int INNER_T=0;INNER_T<q_count;INNER_T++)
//            {
//                double q = random.nextDouble();
//                int query_rank1 = (int) Math.floor(q * (queryN-1)+1),query_rank2 = (int) Math.ceil(q * (queryN-1)+1);
////                System.out.println("\n\t\t\t\t\t\tq rank1,2:"+q+"\t"+query_rank1+","+query_rank2);
//                double[][] deterministic_results=new double[BATCH][];
//                IntArrayList remaining = new IntArrayList();
//                for(int i=0;i<BATCH;i++){
//                    deterministic_results[i]=(lazyWorkers.get(i).findResultRange(query_rank1,query_rank2));
//                    avg_iteration[i]+=ratioPerQuery;
//                    remaining.add(i);
//                }
//                int MMP=0;
//                while(true) {
////                    System.out.println("\t\t"+remaining);
//                    for(int i:new IntArrayList(remaining))
//                        if(deterministic_results[i][0]<deterministic_results[i][1]&&deterministic_results[i].length!=3);
//                    else remaining.removeIf((x)->x==i);
//                    if(remaining.isEmpty())break;
//                    if(++MMP>10){
//                        int id=remaining.getInt(0);
//                        System.out.println("\t\t\t\t\t\titerate fail. too many iter \titer:"+MMP+"\t\tq:"+q+"\t\tMem:\t"+MM.getInt(id)+"\tvalLR:"+ Arrays.toString(deterministic_results[id]));
//                        break;
//                    }
////                    System.out.println("\n\tnew iter.??\t"+deterministic_result[0]+"..."+deterministic_result[1]+"\t\tFilterSize:\t"+last_n);
//
//                    ObjectArrayList<MRLSketchLazy> cntWorkers=new ObjectArrayList<>();
//                    DoubleArrayList valLs=new DoubleArrayList(),valRs=new DoubleArrayList();
//                    int[] CountOfLessThanValLs=new int[remaining.size()],CountOfValLs=new int[remaining.size()],CountOfValRs=new int[remaining.size()];
//                    for(int i:remaining){
//                        avg_iteration[i]+=ratioPerQuery;
//                        cntWorkers.add(new MRLSketchLazy(MM.getInt(i)));
//                        valLs.add(deterministic_results[i][0]);
//                        valRs.add(deterministic_results[i][1]);
//                    }
//                    for(int p=L;p<R;p++) {
//                        double v = a[p];
//                        for(int i=0;i<remaining.size();i++){
//                            double valL=valLs.getDouble(i),valR=valRs.getDouble(i);
//                            if (v > valL && v < valR) cntWorkers.get(i).update(dataToLong(v));
//                            else if (v < valL) CountOfLessThanValLs[i]++;
//                            else if(v==valL)CountOfValLs[i]++;
//                            else if(v==valR)CountOfValRs[i]++;
//                        }
//                    }
//
//
//                    for(int i=0;i<remaining.size();i++){
//                        double valL=valLs.getDouble(i),valR=valRs.getDouble(i);
//                        int CountOfLessThanValL=CountOfLessThanValLs[i],CountOfValL=CountOfValLs[i],CountOfValR=CountOfValRs[i];
//                        MRLSketchLazy cntWorker=cntWorkers.get(i);
//                        int cntRank1=query_rank1-CountOfLessThanValL;
//                        int cntRank2=query_rank2-CountOfLessThanValL;
//
//                        if(cntRank1<=0||cntRank2>CountOfValL+CountOfValR+cntWorker.getN()){ // iteration failed.
//                            System.out.println("\t\t\t\t\t\titerate fail. error filter");
//                            avg_iteration[remaining.getInt(i)]+=10;
//                            remaining.removeInt(i);
//                        }else
//                        deterministic_results[remaining.getInt(i)] = cntWorker.getFilter(CountOfValL,CountOfValR,valL,valR,cntRank1,cntRank2);
//                    }
//                }
//            }
//        }
////        System.out.println("\t\t"+avg_iteration);
//        return avg_iteration;
//    }


    public static boolean noNeed(int i,int j){
//        return /*(i < mNum / 4 && j >= nNum / 4 * 3 + i) || */(i >= mNum / 4 * 3 && j < (i - mNum / 4 * 3 + 1));  //n<=1e8
//        return (i < mNum / 4 && j > i); // 1e8<=n<=1e9
        boolean noNeed=true;
        if(i==20&&j==59)noNeed=false;
        if(i==24&&j==63)noNeed=false;
        if(i==24&&j==64)noNeed=false;
        if(i==25&&j==64)noNeed=false;
        return noNeed;
    }
    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForExactPrintTableGK main;
        main = new MainForExactPrintTableGK();


        IntArrayList NN = new IntArrayList();
        for (double i = 0, n = nST, rela = Math.exp(Math.log(nED / nST) / (nNum - 1)); i < nNum; i++, n *= rela)
            NN.add((int) Math.round(n));
        System.out.println(NN);
        IntArrayList MM = new IntArrayList();
        for (double i = 0, m = mST, rela = Math.exp(Math.log(mED / mST) / (mNum - 1)); i < mNum; i++, m *= rela)
            MM.add((int) Math.round(m));
        System.out.println(MM);


        double[][] MNPass = new double[(int) mNum][(int) nNum];
        main.prepareA(0);
        IntArrayList[] tmm= new IntArrayList[(int) nNum];
        for (int j = 0; j < nNum; j++)tmm[j]=new IntArrayList();

        for (int i = 0; i < mNum; i++)
            for (int j = 0; j < nNum; j++)
                if (noNeed(i,j))
                    MNPass[i][j] = 0;
                else {
                    tmm[j].add(MM.getInt(i));
                }
        for (int j = 0; j < nNum; j++) {
            System.out.println("\t\ttesting line id:\t"+j);
            double[] res = main.testErrorBatch(NN.getInt(j), tmm[j]);
            int tp=0;
            for (int i = 0; i < mNum; i++)
                if (noNeed(i,j));
            else MNPass[i][j]=res[tp++];
        }

//        for (int i = 0; i < mNum; i++)
//            for (int j = 0; j < nNum; j++)
//                if ((i < mNum / 4 && j >= nNum / 4 * 3 + i) || (i >= mNum / 4 * 3 && j < (i - mNum / 4 * 3 + 1)))
//                    MNPass[i][j] = 0;
//                else {
//                    double tmp = main.testError(NN.getInt(j), MM.getInt(i));
//                    MNPass[i][j] = tmp;
//                    System.out.println("\t\t" + MM.getInt(i) + "\t" + NN.getInt(j) + "\t\t" + tmp);
//                }

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

        System.out.println("\t\t\tALL_TIME:" + (new Date().getTime() - START_T));
    }
}