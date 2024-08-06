package rbbeast.beauti;





import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.inference.*;
import beast.base.inference.operator.kernel.BactrianDeltaExchangeOperator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beastfx.app.inputeditor.*;
import beastfx.app.inputeditor.InputEditor.ExpandOption;
import beastfx.app.util.FXUtils;
import javafx.scene.control.CheckBox;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import rbbeast.evolution.sitemodel.FreeRateModel;

import java.util.ArrayList;
import java.util.List;

public class FreeRateModelInputEditor extends SiteModelInputEditor {

	@Override
    public Class<?> type() {
        return FreeRateModel.class;
    }
    
	public FreeRateModelInputEditor(BeautiDoc doc) {
		super(doc);
	}

	@Override
	public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {
		super.init(input, plugin, itemNr, bExpandOption, bAddButtons);
//		InputEditor substModelInputEditor = doc.inputEditorFactory.createInputEditor(((FreeRateModel) input.get()).substModelInput, m_beastObject, doc);
//		
//		pane.getChildren().add(substModelInputEditor);
//
//		validateInput();
	}
	
}
