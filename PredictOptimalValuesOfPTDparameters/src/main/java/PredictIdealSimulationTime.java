import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Scanner;

public class PredictIdealSimulationTime {
	List<IdealTime> idealTimes = null;

	public PredictIdealSimulationTime(String csvPathName) throws FileNotFoundException {
		idealTimes = new ArrayList<IdealTime>();

		if (!new File(csvPathName).exists()) {
			throw new FileNotFoundException(csvPathName + " not exists!");
		}

		Scanner scanner = new Scanner(new FileInputStream(csvPathName));

		scanner.nextLine();

		String line;
		String[] columns;
		while (scanner.hasNextLine()) {
			line = scanner.nextLine();
			columns = line.split(",");

			idealTimes.add(new IdealTime(Integer.parseInt(columns[0]), Double.parseDouble(columns[1])));
			Collections.sort(idealTimes, (p1, p2) -> p1.nProc < p2.nProc ? -1 : 1);
		}
		
		scanner.close();
	}
	
	synchronized double getPredictedTime(double nProc)
	{
		double time = Double.MAX_VALUE;
		
		double x0 , y0;
		double x1 , y1;
		x0 = x1 = y0 = y1 = 0;
		
		IdealTime id1 = null;
		IdealTime id2 = null;
		
		if(idealTimes.get(0).nProc > nProc)
		{
			id1 = idealTimes.get(0);
			id2 = idealTimes.get(1);
		}
		else if(idealTimes.get(idealTimes.size() - 1).nProc < nProc)
		{
			id1 = idealTimes.get(idealTimes.size() - 2);
			id2 = idealTimes.get(idealTimes.size() - 1);
		}
		else
			for (int i = 1; i < idealTimes.size(); i++) 
			{
				id1 = idealTimes.get(i);
				
				if(id1.nProc > nProc)
				{
					id2 = idealTimes.get(i - 1);
					break;
				}
			}

		time = ( id1.time * (id2.nProc - nProc) + id2.time * (nProc - id1.nProc)) / (id2.nProc - id1.nProc );
		
		return time;
	}
	
	@Override
	public String toString() {
		return idealTimes.toString();
	}
}

class IdealTime {
	int nProc;
	double time;

	public IdealTime(int nProc, double time) {
		super();
		this.nProc = nProc;
		this.time = time;
	}

	@Override
	public String toString() {
		return " \n [ " + nProc + ", " + time + " ]";
	}
}
