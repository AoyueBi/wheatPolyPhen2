/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package format.dna;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashByteByteMaps;

/**
 * Utilities related to DNA sequence
 * @author Fei Lu
 */
public class DNAUtils {
    /**The byte value of 4 DNA bases, A, C, G, T*/
    private static final byte[] baseByte = {65, 67, 71, 84};
    /**The byte value of 4 DNA bases and any base, A, C, G, N, T*/
    private static final byte[] baseByteWithN = {65, 67, 71, 78, 84};
    
    private static final byte[] compleBaseByteWithN = {84, 71, 67, 78, 65};
    
    private static final byte[] lowerBaseByte = {97, 99, 103, 116};
    
    /**
     * Return byte value hash map pointing to complementary bases
     * @return 
     */
    public static HashByteByteMap getBaseCompleByteMap () {
        return HashByteByteMaps.getDefaultFactory().newImmutableMap(baseByteWithN, compleBaseByteWithN);
    }
    
    /**
     * Return a byte value base array of A, C, G, T
     * @return 
     */
    public static byte[] getBaseByteArray () {
        return baseByte;
    }
    
    /**
     * Return a byte value base array of A, C, G, N, T
     * @return 
     */
    public static byte[] getBaseWithNByteArray () {
        return baseByteWithN;
    }
}
