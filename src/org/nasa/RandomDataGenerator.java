package org.nasa;

import java.io.BufferedWriter;
import java.io.FileWriter;

public class RandomDataGenerator {

	public static void main(String args[]) throws Exception{
		 FileWriter fstream = new FileWriter("Output1");
		 BufferedWriter out = new BufferedWriter(fstream);
		 int Random = (int)Math.floor(Math.random()*9+1);
		 int j=0;
		while(j<100){ 
			 String result = "";
			 int i=0;
		  while(i<200){
			 	  result = result +(int)Math.floor(Math.random()*10+1) + " " ;
			      i++;		 
		  }
		  System.out.println(result);
		  result = result.substring(0, result.length()-1);
		  out.write(result);
		  out.write("\n");
		  j++ ;
			
		}
		out.close();
		}
	}

