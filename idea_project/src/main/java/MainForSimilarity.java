import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.longs.LongIntPair;
import it.unimi.dsi.fastutil.objects.Object2IntArrayMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Comparator;
import java.util.Date;
import java.util.Objects;


public class MainForSimilarity {
    static int SketchSizeByte = 1024 * 4, SimQuantileNum = 99;

    static DecimalFormat dfPr = new DecimalFormat("#.##E0");
    static DecimalFormat dfPass = new DecimalFormat("0.00");
    static DecimalFormat dfID = new DecimalFormat("00");

    static Object2IntArrayMap<String> fileName2Id = new Object2IntArrayMap<>();
    static ObjectArrayList<String> fileName = new ObjectArrayList<>();
    ObjectArrayList<DoubleArrayList> timeseriesData = new ObjectArrayList<>();
    ObjectArrayList<DoubleArrayList> timeseriesDataSorted = new ObjectArrayList<>();
    ObjectArrayList<KLLSketchLazyExactPriori> timeseriesSketch = new ObjectArrayList<>();
    static double[][] timeseriesDiffKLL, timeseriesDiffData,timeseriesExtJaccardKLL,timeseriesExtJaccardData;

    private String FileToName(String full){
        if(full.indexOf('.')==-1)return "?";
        String a= full.substring(0,full.indexOf('_'));
        String b= full.substring(full.indexOf('_')+1,full.indexOf('.'));
        String c = dfID.format(Integer.parseInt(b));
        return a+'_'+c;
    }

    public void prepareA(String folder) throws IOException {
        ObjectArrayList<File> files = new ObjectArrayList<>(Objects.requireNonNull(new File(folder).listFiles()));
        files.sort((x, y) -> (String.CASE_INSENSITIVE_ORDER.compare(FileToName(x.getName()), FileToName(y.getName()))));
        System.out.println("\t\t\t\tnumber of files:\t" + files.size());
        int fileId = 0;
        for (File file : files) if(!file.isDirectory()){
            String name = file.getName().substring(0, file.getName().lastIndexOf('.'));
            fileName.add(name);
            fileName2Id.put(name, fileId++);
            BufferedReader reader;
            reader = new BufferedReader(new FileReader(file));
            reader.readLine(); // ignore first line.
            String line;
            DoubleArrayList data = new DoubleArrayList();
            while ((line = reader.readLine()) != null)
//                data.add(fileName.size());
                data.add(Double.parseDouble(line));
            timeseriesData.add(data);
            System.out.println("\t\tfile_name:\t" + file + "\t\tlength:\t" + data.size());
            DoubleArrayList dataSorted = new DoubleArrayList(data);
            dataSorted.sort(Double::compare);
            timeseriesDataSorted.add(dataSorted);
        }

        timeseriesDiffKLL = new double[fileName2Id.size()][fileName2Id.size()];
        timeseriesDiffData = new double[fileName2Id.size()][fileName2Id.size()];
        timeseriesExtJaccardKLL = new double[fileName2Id.size()][fileName2Id.size()];
        timeseriesExtJaccardData = new double[fileName2Id.size()][fileName2Id.size()];
    }

    private long dataToLong(double data) {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }

    private double longToResult(long result) {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
    }

    public int getValueActualRank(DoubleArrayList sorted, double v) { // number of elements <= v
        if (v <= sorted.getDouble(0)) return 0;
        int L = 0, R = sorted.size() - 1;
        while (L < R) {
            int mid = (L + R + 1) >>> 1;
            if (v < sorted.getDouble(mid)) R = mid - 1;
            else L = mid;
        }
        return L + 1;
    }

    public void prepareSketch() {
        for (int series = 0; series < timeseriesData.size(); series++) {
            KLLSketchLazyExactPriori sketch = new KLLSketchLazyExactPriori(SketchSizeByte);
            for (double value : timeseriesData.get(series))
                sketch.update(dataToLong(value));
            timeseriesSketch.add(sketch);
        }
    }

    public void computeDisSimBySketch() {
        for (int seriesA = 0; seriesA < timeseriesData.size(); seriesA++)
            for (int seriesB = 0; seriesB < timeseriesData.size(); seriesB++) /*if(seriesA!=seriesB)*/ {
                KLLSketchLazyExactPriori sketchA = timeseriesSketch.get(seriesA), sketchB = timeseriesSketch.get(seriesB);
//                if(seriesA!=15+9 || (seriesB != 0&&seriesB!=15))continue;
//                sketchA.showNum();
//                sketchB.showNum();
                double diff = 0;
//                for (int q_th = 1; q_th <= SimQuantileNum; q_th += 1) {
//                    double q = 1.0 * q_th / (SimQuantileNum + 1);
//                    long longValA = sketchA.findMinValueWithRank((long) (q * sketchA.getN()));
//                    long rankB = sketchB.getApproxRank(longValA);
//                    double q_B = 1.0 * rankB / sketchB.getN();
////                    System.out.println("\t\t"+Math.abs(q - q_B)+"\t\t"+longToResult(longValA)+"\t"+q+","+q_B);
//                    diff += Math.abs(q - q_B) / SimQuantileNum;
//                }
//                timeseriesDiffKLL[seriesA][seriesB] = diff;

                diff = 0;
                LongArrayList rankAs=new LongArrayList();
                for (int q_th = 1; q_th <= SimQuantileNum; q_th += 1) rankAs.add((long) (1.0 * q_th / (SimQuantileNum + 1) * sketchA.getN()));
                LongArrayList longValAs = sketchA.findMinValueWithRankMulti(rankAs,sketchA.getMin(),sketchB.getMax());
                for (int q_th = 1; q_th <= SimQuantileNum; q_th += 1) {
                    double q = 1.0 * q_th / (SimQuantileNum + 1);
                    long longValA = longValAs.getLong(q_th-1);
                    long rankB = sketchB.getApproxRank(longValA);
                    double q_B = 1.0 * rankB / sketchB.getN();
                    diff += Math.abs(q - q_B) / SimQuantileNum;
                }
//                System.out.println("\t\t\t\t"+timeseriesDiffKLL[seriesA][seriesB] +"\t\t"+diff);
                timeseriesDiffKLL[seriesA][seriesB] = diff;
            }
    }

    public void computeSetExtJaccardBySketch() {
        for (int seriesA = 0; seriesA < timeseriesData.size(); seriesA++)
            for (int seriesB = 0; seriesB < timeseriesData.size(); seriesB++) /*if(seriesA!=seriesB)*/ {
                KLLSketchLazyExactPriori sketchA = timeseriesSketch.get(seriesA), sketchB = timeseriesSketch.get(seriesB);
                long[] numA = sketchA.getNum(),numB= sketchB.getNum();
                int[] lpA=sketchA.getLevelPos(),lpB=sketchB.getLevelPos();
                ObjectArrayList<LongIntPair> itemA = new ObjectArrayList<>(),itemB = new ObjectArrayList<>(), itemUnion = new ObjectArrayList<>(), itemInter = new ObjectArrayList<>();
                for(int lv=0;lv<sketchA.getCntLevel();lv++)for(int i=lpA[lv];i<lpA[lv+1];i++)itemA.add(LongIntPair.of(numA[i],1<<lv));
                for(int lv=0;lv<sketchB.getCntLevel();lv++)for(int i=lpB[lv];i<lpB[lv+1];i++)itemB.add(LongIntPair.of(numB[i],1<<lv));
                itemA.sort(Comparator.comparingLong(LongIntPair::firstLong));
                itemB.sort(Comparator.comparingLong(LongIntPair::firstLong));
                ObjectArrayList<LongIntPair> itemAB = new ObjectArrayList<>();
                itemAB.addAll(itemA);
                itemAB.addAll(itemB);
                itemAB.sort(Comparator.comparingLong(LongIntPair::firstLong));
                int extUnion = 0, extInter=0,posA=0,posB=0;
                for(int i=0;i<itemAB.size();i++)if(i==0||itemAB.get(i).firstLong()!=itemAB.get(i-1).firstLong()){
                    long cntV = itemAB.get(i).firstLong();
                    int extA =0,extB=0;
                    while(posA<itemA.size()&&itemA.get(posA).firstLong()==cntV)extA+=itemA.get(posA++).secondInt();
                    while(posB<itemB.size()&&itemB.get(posB).firstLong()==cntV)extB+=itemB.get(posB++).secondInt();
                    extUnion+=Math.max(extA,extB);
                    extInter+=Math.min(extA,extB);
                }
                timeseriesExtJaccardKLL[seriesA][seriesB] = 1.0*extInter/extUnion;
            }
    }

    public void computeDisSimByData() {
        for (int seriesA = 0; seriesA < timeseriesData.size(); seriesA++)
            for (int seriesB = 0; seriesB < timeseriesData.size(); seriesB++) /*if(seriesA!=seriesB)*/ {
//                DoubleArrayList dataSortedA = timeseriesDataSorted.get(seriesA), dataSortedB = timeseriesDataSorted.get(seriesB);
                DoubleArrayList dataSortedA = new DoubleArrayList(timeseriesData.get(seriesA)), dataSortedB = new DoubleArrayList(timeseriesData.get(seriesB));
                dataSortedA.sort(Double::compare);
                dataSortedB.sort(Double::compare);

//                if(seriesA!=15+9 || (seriesB != 0&&seriesB!=15))continue;
//                System.out.println("\t\t"+dataSortedA);
//                System.out.println("\t\t"+dataSortedB);
                double diff = 0;
                int cntRankB = 0;
                for (int q_th = 1; q_th <= SimQuantileNum; q_th += 1) {
                    double q = 1.0 * q_th / (SimQuantileNum + 1);
                    double valA = dataSortedA.getDouble((int) (q * dataSortedA.size()));
                    while (cntRankB < dataSortedB.size() && valA >= dataSortedB.getDouble(cntRankB)) cntRankB++;
                    double q_B = 1.0 * cntRankB / dataSortedB.size();
                    diff += Math.abs(q - q_B) / SimQuantileNum;
                }
                timeseriesDiffData[seriesA][seriesB] = diff;
            }
    }



    public void computeSetExtJaccardByData() {
        for (int seriesA = 0; seriesA < timeseriesData.size(); seriesA++)
            for (int seriesB = 0; seriesB < timeseriesData.size(); seriesB++) /*if(seriesA!=seriesB)*/ {
//                DoubleArrayList dataSortedA = timeseriesDataSorted.get(seriesA), dataSortedB = timeseriesDataSorted.get(seriesB);
                DoubleArrayList dataSortedA = new DoubleArrayList(timeseriesData.get(seriesA)), dataSortedB = new DoubleArrayList(timeseriesData.get(seriesB));
                dataSortedA.sort(Double::compare);
                dataSortedB.sort(Double::compare);

                DoubleArrayList dataSortedAB = new DoubleArrayList();
                for (int pA = 0, pB = 0; pA < dataSortedA.size() || pB < dataSortedB.size(); )
                    dataSortedAB.add((pB == dataSortedB.size() || (pA < dataSortedA.size() && dataSortedA.getDouble(pA) <= dataSortedB.getDouble(pB))) ? dataSortedA.getDouble(pA++) : dataSortedB.getDouble(pB++));

                int extUnion = 0, extInter = 0, posA = 0, posB = 0;
                for (int i = 0; i < dataSortedAB.size(); i++)
                    if (i == 0 || dataSortedAB.getDouble(i) != dataSortedAB.getDouble(i - 1)) {
                        double cntV = dataSortedAB.getDouble(i);
                        int extA = 0, extB = 0;
                        while (posA < dataSortedA.size() && dataSortedA.getDouble(posA) == cntV) {extA++;posA++;}
                        while (posB < dataSortedB.size() && dataSortedB.getDouble(posB) == cntV) {extB++;posB++;}
                        extUnion += Math.max(extA, extB);
                        extInter += Math.min(extA, extB);
                    }
                timeseriesExtJaccardData[seriesA][seriesB] = 1.0 * extInter / extUnion;
            }
    }

    public void printDistributionDiffMatrix(double[][] resultKLL,double[][] resultData){
        System.out.println("\n\t[Distribution Diff (Quantile number "+SimQuantileNum+")] by sketch (approx quantile and rank):");
        System.out.print("\t");
        for (int i = 0; i < fileName.size(); i++) System.out.print("\t" + fileName.get(i));
        System.out.println();

        for (int i = 0; i < fileName.size(); i++) {
            System.out.print("\t" + fileName.get(i) + "\t" + i);
            for (int j = 0; j < fileName.size(); j++)
                System.out.print("\t" + resultKLL[i][j]);
            System.out.println();
        }
        System.out.println("\n\t[Distribution Diff (Quantile number \"+SimQuantileNum+\")] by data:");
        System.out.print("\t");
        for (int i = 0; i < fileName.size(); i++) System.out.print("\t" + fileName.get(i));
        System.out.println();
        for (int i = 0; i < fileName.size(); i++) {
            System.out.print("\t" + fileName.get(i) + "\t" + i);
            for (int j = 0; j < fileName.size(); j++)
                System.out.print("\t" + resultData[i][j]);
            System.out.println();
        }
    }

    public void printJACCARDMatrix(double[][] resultKLL,double[][] resultData){
        System.out.println("\n\t[Extended Jaccard] by sketch:");
        System.out.print("\t");
        for (int i = 0; i < fileName.size(); i++) System.out.print("\t" + fileName.get(i));
        System.out.println();
        for (int i = 0; i < fileName.size(); i++) {
            System.out.print("\t" + fileName.get(i) + "\t" + i);
            for (int j = 0; j < fileName.size(); j++)
                System.out.print("\t" + resultKLL[i][j]);
            System.out.println();
        }
        System.out.println("\n\t[Extended Jaccard] by data:");
        System.out.print("\t");
        for (int i = 0; i < fileName.size(); i++) System.out.print("\t" + fileName.get(i));
        System.out.println();
        for (int i = 0; i < fileName.size(); i++) {
            System.out.print("\t" + fileName.get(i) + "\t" + i);
            for (int j = 0; j < fileName.size(); j++)
                System.out.print("\t" + resultData[i][j]);
            System.out.println();
        }
    }


    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForSimilarity main = new MainForSimilarity();

        main.prepareA("./dataset");
        main.prepareSketch();

        long TIME_DIFF_SKETCH = -new Date().getTime();
        main.computeDisSimBySketch();
        TIME_DIFF_SKETCH+=new Date().getTime();

        long TIME_JACCARD_SKETCH = -new Date().getTime();
        main.computeSetExtJaccardBySketch();
        TIME_JACCARD_SKETCH+=new Date().getTime();

        long TIME_DIFF_DATA = -new Date().getTime();
        main.computeDisSimByData();
        TIME_DIFF_DATA+=new Date().getTime();
        long TIME_JACCARD_DATA = -new Date().getTime();
        main.computeSetExtJaccardByData();
        TIME_JACCARD_DATA+=new Date().getTime();

        main.printDistributionDiffMatrix(timeseriesDiffKLL,timeseriesDiffData);

        main.printJACCARDMatrix(timeseriesExtJaccardKLL,timeseriesExtJaccardData);

        System.out.println("\tDIS_DIFF_TIME:\t\tsketch:\t" + TIME_DIFF_SKETCH);
        System.out.println("\tDIS_DIFF_TIME:\t\tdata:\t" + TIME_DIFF_DATA);
        System.out.println("\tDIS_JACCARD_TIME:\t\tsketch:\t" + TIME_JACCARD_SKETCH);
        System.out.println("\tDIS_JACCARD_TIME:\t\tdata:\t" + TIME_JACCARD_DATA);

        System.out.println("\t\t\tALL_TIME:" + (new Date().getTime() - START_T));
    }
}
