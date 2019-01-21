package explore;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

public class RelativizePaths
{
	public static void main(String[] args)
	{
		Path targetPath = Paths.get("/var/data/stuff/xyz.dat");

		ArrayList< Path > paths = new ArrayList<>(  );
		paths.add( Paths.get("/var/data") );
		paths.add( Paths.get("/var/data/other") );
		paths.add( Paths.get("/var/data/other/file.dat") );

		for ( Path path : paths )
		{
			System.out.println( path.relativize( targetPath ) );
		}
	}
}
