package org.nasa;
import java.util.*;

public class test {
    static Map <String,Integer> dataMap = new HashMap<String,Integer>();
	static{
		
		char a = 'a' ;
		
		for(int i=0;i<26;i++){
			String temp =  "" +((char)(a+i));
			System.out.println(temp);
			dataMap.put(temp, i);
		}
		
	}
	
	public static int getInteger(String obs){
		return dataMap.get(obs);
	}
}
