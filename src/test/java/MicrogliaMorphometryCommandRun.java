import de.embl.cba.morphometry.commands.MicrogliaMorphometryCommand;
import de.embl.cba.morphometry.commands.TranslocationCommand;
import ij.IJ;
import ij.plugin.frame.RoiManager;
import net.imagej.ImageJ;

public class MicrogliaMorphometryCommandRun
{
	public static void main( String[] args )
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		ij.command().run( MicrogliaMorphometryCommand.class, true );
	}
}
