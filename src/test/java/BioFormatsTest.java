import loci.formats.FilePattern;

import java.io.File;

import automic.table.TableModel;
import automic.utils.imagefiles.PatternOpenerWithBioformats;
import ij.ImageJ;
import ij.ImagePlus;

public abstract class BioFormatsTest {


	static final String tablePath="/Volumes/almfscreen/kletter/20180904/formattedLog_0906-1_Sorted_Cropped.txt";
	static final String highZoomPatternTag="High.Zoom";

	public static void main(String[] args)throws Exception{

		//specify folder with your ImageJ plugins. This allows IJ.run("Bio-Formats",...) to work
		System.getProperties().setProperty("plugins.dir", "/Applications/Fiji.app/plugins");

		//test imageJ instance
		new ImageJ();

		//create table model from file
		TableModel tbl=new TableModel(new File(tablePath));

		//get number of datasets in the table
		int nDatasets=tbl.getRowCount();

		//loop through datasets
		for (int iDataset=0; iDataset<nDatasets; iDataset++) {
			if(!tbl.getBooleanValue(iDataset, "Success"))//no high zoom image acquired
				continue;

			//get folder and filename components of the pattern from the table
			String patternPath=tbl.getFileParentPathString(iDataset, highZoomPatternTag,"IMGR");
			String patternName=tbl.getFileName(iDataset, highZoomPatternTag,"IMGR");

			ImagePlus highZoomImage=PatternOpenerWithBioformats.openImage(new File(patternPath,patternName));

			//show image for 5 seconds
			highZoomImage.show();
			Thread.sleep(5000);
			highZoomImage.hide();

		}
	}
}
