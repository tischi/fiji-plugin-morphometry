package de.embl.cba.morphometry;

import ij.IJ;

import java.io.File;

public class Logger
{
	public static void debug( String s )
	{
		//IJ.log( s );
	}

	public static void log( String message )
	{
		IJ.log( message );

		if ( Utils.logFilePath != null )
		{
			File logFile = new File( Utils.logFilePath );

			if ( ! logFile.exists() )
			{
				Utils.createLogFile();
			}
			else
			{
				Utils.writeToLogFile( message + "\n" );
			}
		}

	}

	public static void error( String s )
	{
		IJ.showMessage( s );
	}
}
