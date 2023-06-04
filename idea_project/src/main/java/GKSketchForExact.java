import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.Serializable;
import java.util.Arrays;

/** Util for UDAFPercentile */
public class GKSketchForExact {
  private final double rankAccuracy;
  private ObjectArrayList<Tuple> entries;
  private int maxEntryNum,entryNum;

  private final double[] buffer;
  private int maxBufferNum, bufferNum;
  private long N;

  public GKSketchForExact(int maxMemoryByte) {
    maxBufferNum=maxMemoryByte/8/8;
    maxEntryNum=(maxMemoryByte/8-maxBufferNum)/3;
    this.rankAccuracy = 1.0/(maxBufferNum-1);
    this.entries = new ObjectArrayList<>(maxBufferNum);
    this.buffer = new double[maxBufferNum];
    this.bufferNum = 0;
    this.N = 0;
  }

  public void insert(double value) {
    buffer[bufferNum++] = value;
    if (bufferNum == maxBufferNum) {
      compress();
    }
  }

  public boolean isEmpty() {
    return entries.isEmpty() && bufferNum == 0;
  }

  public double query(double phi) {
    if (isEmpty()) {
      throw new ArithmeticException();
    }
    compressIfNecessary();

    long maxGD=0;
    for(Tuple entry:entries)maxGD=Math.max(maxGD,entry.g+entry.delta);
    long targetErr = maxGD/2;
    final long rank = (long) (phi * (N - 1)) + 1;
    long gSum = 0,maxRank;
    for (Tuple entry : entries) {
      gSum += entry.g;
      maxRank = gSum + entry.delta;
      if (maxRank - targetErr <= rank && rank <= gSum + targetErr) {
        return entry.v;
      }
    }
    return entries.get(entries.size()-1).v;
  }

  private void compressIfNecessary() {
    if (bufferNum > 0) {
      compress();
    }
  }

  private void compress() {
//    N += bufferNum;
    Arrays.sort(buffer,0, bufferNum);

    final ObjectArrayList<Tuple> mergedEntries =
        new ObjectArrayList<>(entryNum + bufferNum);

    int pe=entryNum-1,pb=bufferNum-1;
    while(pe>=0||pb>=0){
      if(pb<0||(pe>=0&&entries.get(pe).v>buffer[pb])){
        mergedEntries.add(entries.get(pe));
        pe--;
      }else{
        N++;
        long delta = (pb==0&&pe<0)||(pe==entryNum-1&&pb==bufferNum-1) ? 0:(long)Math.floor(2*rankAccuracy*N);
        mergedEntries.add(new Tuple(buffer[pb],1,delta));
        pb--;
      }
    }
    entries.clear();
    for(int i=mergedEntries.size()-1;i>=0;i--)entries.add(mergedEntries.get(i));
    entryNum=entries.size();
    bufferNum = 0;

    if(entryNum<=maxEntryNum)return;

    final double removalThreshold = (2 * rankAccuracy *N);
//    System.out.print("\tbeforeCompress\tval(g+delta):\t");
//    for(Tuple entry:entries)System.out.print("\t"+entry.v+"("+(entry.g+entry.delta)+")");
//    System.out.println();

    Tuple last = entries.get(entryNum-1),now;
    final ObjectArrayList<Tuple> revEntries =
        new ObjectArrayList<>(entryNum + bufferNum);
    for(int i=entryNum-2;i>=1;i--){
      now=entries.get(i);
      if(now.g+last.g+last.delta<removalThreshold)
        last.g+=now.g;
      else{
        revEntries.add(last);
        last=now;
      }
    }
    revEntries.add(last);
    if(entries.get(0).v<=last.v&&entryNum>1)
      revEntries.add(entries.get(0));

    entries.clear();
    for(int i=revEntries.size()-1;i>=0;i--)entries.add(revEntries.get(i));
    entryNum=entries.size();

//    if(1.0*entryNum/maxEntryNum>=1.0)
//      show();
//      System.out.println("N:"+N+" entry:"+entryNum+"\t"+1.0*entryNum/maxEntryNum);
//    System.out.print("N:"+N+" entry:"+entryNum+"\t"+1.0*entryNum/maxEntryNum+"\tThreshold:\t"+removalThreshold+"\t\tval(g+delta):\t");
//    for(Tuple entry:entries)System.out.print("\t"+entry.v+"("+(entry.g+entry.delta)+")");
//    System.out.println();
  }

  public void show(){
      System.out.print("N:"+N+" entry:"+entryNum+"\t"+1.0*entryNum/maxEntryNum+"\t\tval(g+delta):\t");
    for(Tuple entry:entries)System.out.print("\t"+entry.v+"("+(entry.g+entry.delta)+")");
    System.out.println();
  }


  public double[] getFilterL(long CountOfValL,long CountOfValR,double valL,double valR,long K){
    if(K<=CountOfValL)return new double[]{valL,-233.0};
    if(K>CountOfValL+N)return new double[]{valR,-233.0};
    K-=CountOfValL;
    double lb=entries.get(0).v;//=minValueWithRank(K-1);

    long gSum = 0,maxRank;
    for (Tuple entry : entries) {
      gSum += entry.g;
      maxRank = gSum + entry.delta;
      if(maxRank>=K){break;}
      lb=entry.v;
    }
    return new double[]{lb};
  }

  public double[] getFilterR(long CountOfValL,long CountOfValR,double valL,double valR,long K){
    if(K<=CountOfValL)return new double[]{valL,-233.0};
    if(K>CountOfValL+N)return new double[]{valR,-233.0};
    K-=CountOfValL;
    double ub=entries.get(entries.size()-1).v;//maxValueWithRank(K);
    long gSum = 0,maxRank;
    for (Tuple entry : entries) {
      gSum += entry.g;
      maxRank = gSum + entry.delta;
      if(gSum>=K){ub=entry.v;break;}
    }
    return new double[]{ub};
  }


  public double[] getFilter(long CountOfValL,long CountOfValR,double valL,double valR,long K1,long K2){
    compressIfNecessary();
    double[] filterL =getFilterL(CountOfValL,CountOfValR,valL,valR,K1);
    double[] filterR =getFilterR(CountOfValL,CountOfValR,valL,valR,K2);
//    System.out.println("\t\t\t\tvalL,R:\t"+filterL[0]+"..."+filterR[0]);
    if(filterL.length+filterR.length==4)return new double[]{filterL[0],filterR[0],-233};
    else return new double[]{filterL[0],filterR[0]};
  }


  public long findMaxNumberInRange(double L,double R){ // L,R is value of entry
    compressIfNecessary();
    long gSum = 0,maxRank;
    long LMinRank=1,RMaxRank=N;
    for (Tuple entry : entries) {
      gSum += entry.g;
      maxRank = gSum + entry.delta;
      if(entry.v==L)LMinRank=gSum;
      if(entry.v==R)RMaxRank=maxRank;
    }
    return RMaxRank-LMinRank+1;
  }

  static class Tuple implements Serializable {
    private final double v;
    public long g;
    public long delta;

    private Tuple(double v, long g, long delta) {
      this.v = v;
      this.g = g;
      this.delta = delta;
    }
  }
}

