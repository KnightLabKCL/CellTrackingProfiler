/**  @author caroff  **/
package plugins.caroff.detection;

import icy.util.XMLUtil;

import org.w3c.dom.Element;
import org.w3c.dom.Node;

/** Detection class to allow the Detection of cells with labels in the Track Manager plugin **/
public class Detection extends plugins.nchenouard.spot.Detection
{
	// Addition of attribute
	double label = 0;
	
	// Modification of the function to load the label of cells
	@Override
	public boolean loadFromXML(Node node)
	{
		super.loadFromXML(node);
		Element detectionElement = (Element) node;
		label = XMLUtil.getAttributeDoubleValue( detectionElement , "label" ,  0 );
		return true;
	}
	
	// Modification of the function to save the label of cells
	@Override
	public boolean saveToXML(Node node)
	{
		super.saveToXML(node);
		Element detectionElement = (Element) node;
		XMLUtil.setAttributeDoubleValue(detectionElement, "label", label );
		return true;
	} 
}
