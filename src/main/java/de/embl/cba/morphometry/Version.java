package de.embl.cba.morphometry;


import org.apache.maven.model.Model;
import org.apache.maven.model.io.xpp3.MavenXpp3Reader;
import org.codehaus.plexus.util.xml.pull.XmlPullParserException;

import java.io.FileReader;
import java.io.IOException;

public class Version
{

	public static String getArtifactVersion()
	{
		MavenXpp3Reader reader = new MavenXpp3Reader();
		Model model = getMavenModel( reader );

//		System.out.println(model.getId());
//		System.out.println(model.getGroupId());
//		System.out.println(model.getArtifactId());
//		System.out.println(model.getVersion());
//
		return model.getVersion();
	}

	public static Model getMavenModel( MavenXpp3Reader reader )
	{
		Model model = null;
		try
		{
			model = reader.read( new FileReader( "pom.xml" ) );
		}
		catch ( IOException e )
		{
			e.printStackTrace();
		}
		catch ( XmlPullParserException e )
		{
			e.printStackTrace();
		}
		return model;
	}
}
