/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analysis.maize2k;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import utils.IOUtils;


/**
 *
 * @author feilu
 */
public class Maize2kGo {
    
    public Maize2kGo () {
        //this.processReference();
        //this.seqQualityTestPipe();
        this.setHapScannerPipe();
        //this.test();
        //today is Apr 9  111   222 333
    }
    
    private void setHapScannerPipe () {
        //new HapMap3Processor();
        String parameterFileS = "/Users/Aoyue/NetBeansProjects/HapScanner/HapScanner/src/analysis/Parameters_hapScanner.txt";
        new HapScanner(parameterFileS);
    }
    
    public void test () {
        String infileS = "/Users/feilu/Documents/analysisL/pipelineTest/HapScanner/hapMap3_AGPV4/hmp321_agpv4_chr10.vcf.gz";
        String outfileS = "/Users/feilu/Documents/analysisL/pipelineTest/HapScanner/test.vcf";
        int length = 10000;
        try {
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < length; i++) {
                bw.write(br.readLine());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void seqQualityTestPipe () {
        new FastqQuality ();
    }
    
    private void processReference () {
        new ReferenceProcessor();
    }
    
    public static void main (String[] args) {
        new Maize2kGo();
    }
    
}
