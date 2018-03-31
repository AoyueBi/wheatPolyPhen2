/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maize2k;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.List;
import utils.IOUtils;
import utils.PStringUtils;


/**
 *
 * @author feilu
 */
public class HapMap3Processor {
    
    public HapMap3Processor () {
        //this.mkHapPosAllele();
        //this.mkHapPos();
        //this.test();
    }
    
    public void mkHapPos () {
        String infileDirS = "/Users/feilu/Documents/analysisL/pipelineTest/HapScanner/hapPosAllele/";
        String outfileDirS = "/Users/feilu/Documents/analysisL/pipelineTest/HapScanner/hapPos/";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".posAllele.txt.gz");
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".posAllele.txt.gz", ".pos.txt.gz")).getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = null;
            int cnt = 0;
            try {
                StringBuilder sb = new StringBuilder();
                List<String> l = null;
                temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    sb = new StringBuilder(l.get(0));
                    sb.append("\t").append(l.get(1));
                    bw.write(sb.toString());
                    bw.newLine();
                    if (cnt%100000 == 0) System.out.println("Output " + String.valueOf(cnt) + " SNPs");
                    cnt++;
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(String.valueOf(cnt) + " SNPs output from " + f.getAbsolutePath());
            }
            catch (Exception e) {
                System.out.println(f.getAbsolutePath());
                System.out.println(temp);
                e.printStackTrace();
                System.exit(1);
            }
        });
    }
    
    public void mkHapPosAllele () {
        String infileDirS = "/Users/feilu/Documents/analysisL/pipelineTest/HapScanner/hapMap3_AGPV4/";
        String outfileDirS = "/Users/feilu/Documents/analysisL/pipelineTest/HapScanner/hapPosAllele/";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".posAllele.txt.gz")).getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = null;
            int cnt = 0;
            try {
                StringBuilder sb = new StringBuilder();
                sb.append("Chr\tPos\tRef\tAlt");
                bw.write(sb.toString());
                bw.newLine();
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) continue;
                    temp = temp.substring(0, 40);
                    l = PStringUtils.fastSplit(temp);
                    sb = new StringBuilder(l.get(0));
                    sb.append("\t").append(l.get(1)).append("\t").append(l.get(3)).append("\t").append(l.get(4));
                    bw.write(sb.toString());
                    bw.newLine();
                    if (cnt%100000 == 0) System.out.println("Output " + String.valueOf(cnt) + " SNPs");
                    cnt++;
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(String.valueOf(cnt) + " SNPs output from " + f.getAbsolutePath());
            }
            catch (Exception e) {
                System.out.println(f.getAbsolutePath());
                System.out.println(temp);
                e.printStackTrace();
                System.exit(1);
            }
            
        });
    }
    
    public void test () {
        
    }
}
