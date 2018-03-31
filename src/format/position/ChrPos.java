/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package format.position;

/**
 *
 * @author feilu
 */
public class ChrPos {
    short chr;
    int pos;
    
    public ChrPos (short chr, int pos) {
        this.chr = chr;
        this.pos = pos;
    }
    
    public short getChr () {
        return chr;
    }
    
    public int getPos() {
        return pos;
    }
}
