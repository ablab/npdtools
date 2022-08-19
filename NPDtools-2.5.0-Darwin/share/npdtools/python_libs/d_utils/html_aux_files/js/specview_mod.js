// $LastChangedDate$
// $LastChangedBy$
// $LastChangedRevision$

$(function() {

    // plugin name - specview_dereplicator
    	$.fn.specview_dereplicator = function (opts) {

        var defaults = {
                sequence: null,
                selectedPeak: null,
                spectrumInfo: null, // dict in format {filename, scan, mz, retention, charge}  //TBD: intensity
                compoundInfo: null, // dict in format {name, formula, mass, adduct}
                spectrum: [], // list of all MS2 peaks in format [m/z, intensity]
                annotation: [], // list of annotated peaks in format {peak_idx, theorMass, charge, components}
                mol: null,
                atomComponents: [],
                modification: null, // dict in format {position, componentFormula, componentMass, massShift} or null if no modification

                width: 600,     // width of the ms/ms plot
                height: 300,    // height of the ms/ms plot
                labelWidth: 35,
                labelMargin: 5,
                showViewingOptions: true,
                showMassErrorPlot: false,
                massErrorUnit: massErrorType_Da,
                precision: 3
        };

        var options = $.extend(true, {}, defaults, opts); // this is a deep copy
        var nonSelectedSpecs, selectedSpecs, modifiedSpecs, modSelectedSpecs, invisibleSpecs, fragmentedBondsSpecs;

        return this.each(function() {

            index = index + 1;
            init($(this), options);

        });
    };

    var index = 0;
    var massErrorType_Da = 'Da';
    var massErrorType_ppm = 'ppm';

    var elementIds = {
            massError: "massError",
            currentMass: "currentMass",
            massErrorPlot: "massErrorPlot",
            massErrorPlot_unit: "massErrorUnit",
            massErrorPlot_option: "massErrorPlot_option",
            msmsplot: "msmsplot",
            ms2plot_zoom_out: "ms2plot_zoom_out",
            zoom_x: "zoom_x",
            zoom_y: "zoom_y",
            resetZoom: "resetZoom",
            saveLink: "saveLink",
            update: "update",
            enableTooltip: "enableTooltip",
            msmstooltip: "lorikeet_msmstooltip",
            slider_width: "slider_width",
            slider_width_val: "slider_width_val",
            slider_height: "slider_height",
            slider_height_val: "slider_height_val",
            psm_general_info: "psm_content",
            lorikeet_content: "lorikeet_content",
            viewOptionsDiv: "viewOptionsDiv",
            peaksTable: "peaksTable",
            peakRow: "peak_row",
            peaksTableDiv: "peaksTableDiv",
            seqinfo: "seqinfo",
            chemCanvas: "chemCanvas",
            seqCanvas: "seq_span"
    };

    var molColors = {
        'default': 'black',
        'selected': 'blue',
        'modified': '#ff9933',
        'mod_selected': '#00ccff'
    };

    function getElementId(container, elementId){
        return elementId+"_"+container.data("index");
    }

    function getRadioName(container, name) {
        return name+"_"+container.data("index");
    }

    function getElementSelector(container, elementId) {
        return "#"+getElementId(container, elementId);
    }

    function getMassErrorPrecision(options) {
        if (options.massErrorUnit == massErrorType_ppm) {
            return 1;  // special case
        }
        else {
            return options.precision;
        }
    }

    function sortByKey(array, key) {
        return array.sort(function(a, b) {
            var x = a[key]; var y = b[key];
            return ((x < y) ? -1 : ((x > y) ? 1 : 0));
        });
    }

    function defineBonds() {
        nonSelectedSpecs = new ChemDoodle.structures.VisualSpecifications();
        nonSelectedSpecs.bonds_color = molColors['default'];
        nonSelectedSpecs.atoms_color = molColors['default'];
        selectedSpecs = new ChemDoodle.structures.VisualSpecifications();
        selectedSpecs.bonds_color = molColors['selected'];
        selectedSpecs.atoms_color = molColors['selected'];
        modifiedSpecs = new ChemDoodle.structures.VisualSpecifications();
        modifiedSpecs.bonds_color = molColors['modified'];
        modifiedSpecs.atoms_color = molColors['modified'];
        modSelectedSpecs = new ChemDoodle.structures.VisualSpecifications();
        modSelectedSpecs.bonds_color = molColors['mod_selected'];
        modSelectedSpecs.atoms_color = molColors['mod_selected'];
        invisibleSpecs = new ChemDoodle.structures.VisualSpecifications();
        invisibleSpecs.bonds_color = 'white';
        invisibleSpecs.atoms_color = 'white';
        fragmentedBondsSpecs = new ChemDoodle.structures.VisualSpecifications();
        fragmentedBondsSpecs.bonds_color = 'red';
    }

    function init(parent_container, options) {
        options.annotation = sortByKey(options.annotation, 'peakIdx');  // i.e. sorting by 'mz'
        var peaksGroups = splitPeaks(options.spectrum, options.annotation);
        options.notMatchedPeaks = peaksGroups[0];
        options.matchedPeaks = peaksGroups[1];

        options.massErrorsDa = calcMassErrors(options.matchedPeaks, options.annotation,
            massErrorType_Da);
        options.massErrorsPpm = calcMassErrors(options.matchedPeaks, options.annotation,
            massErrorType_ppm);
        refreshMassErrors(options);

        if (options.spectrumInfo || options.compoundInfo)
            addGeneralInfoPanel(parent_container, options);

        var container = createContainer(parent_container);
        // alert(container.attr('id')+" parent "+container.parent().attr('id'));
        storeContainerData(container, options);
        initContainer(container);
        defineBonds();

        makeViewingOptions(container, options);

        if(options.sequence) {
            showSequenceInfo(container, options);
        }

        createPlot(container, getDatasets(container)); // Initial MS/MS Plot

        setupInteractions(container, options);
        if (options.mol) showNewMol(options.mol, options);
        addMatchedPeaksTable(container);
    }

    function splitPeaks(spectrum, annotation) {
        var allPeakIdx = 0;
        var matchedPeakIdx = 0;
        var curMatchedPeakInfo = matchedPeakIdx < annotation.length ?
            annotation[matchedPeakIdx] : null;
        var notMatchedPeaks = [];
        var matchedPeaks = [];
        while (allPeakIdx < spectrum.length) {
            if (curMatchedPeakInfo && curMatchedPeakInfo.peakIdx == allPeakIdx) {
                matchedPeaks.push(spectrum[allPeakIdx]);
                matchedPeakIdx++;
                curMatchedPeakInfo = matchedPeakIdx < annotation.length ? annotation[matchedPeakIdx] : null;
            } else {
                notMatchedPeaks.push(spectrum[allPeakIdx]);
            }
            allPeakIdx++;
        }
        return [notMatchedPeaks, matchedPeaks];
    }

    function refreshMassErrors(options) {
        if (options.massErrorUnit == massErrorType_ppm) {
            options.massErrors = options.massErrorsPpm;
        }
        else {
            options.massErrors = options.massErrorsDa;
        }
    }

    function calcMassErrors(peaks, annotation, unit) {  // unit is massErrorType_Da or massErrorType_ppm
        var massErrors = [];
        for (var i = 0; i < peaks.length; i++) {
            var peakInfo = annotation[i];
            var peakMass = peaks[i][0] - Ion.MASS_PROTON; // do not need to mult by charge because it is taken into account in theorMass
            var theorMass = peakInfo.theorMass;
            var massError = peakMass - theorMass;
            if (unit == massErrorType_ppm) {
                massError = massError * Math.pow(10, 6) / theorMass;
            }
            massErrors.push(massError);
        }
        return massErrors;
    }

    function storeContainerData(container, options) {
        container.data("index", index);
        container.data("options", options);
        container.data("massErrorChanged", false);
        container.data("massTypeChanged", false);
        container.data("plot", null);           // MS/MS plot
        container.data("zoomRange", null);      // for zooming MS/MS plot
        container.data("previousPoint", null);  // for tooltips
        container.data("massError", 0.0);

        var maxInt = getMaxInt(options);
        var __xrange = getPlotXRange(options);
        var plotOptions =  {
                series: {},
                selection: { mode: "x", color: "#F0E68C" },
                grid: { show: true,
                        hoverable: true,
                        autoHighlight: false,
                        clickable: true,
                        borderWidth: 1,
                        labelMargin: options.labelMargin },
                xaxis: { tickLength: 3, tickColor: "#000",
                         min: __xrange.xmin,
                         max: __xrange.xmax},
                yaxis: { tickLength: 0, tickColor: "#000",
                         max: maxInt*1.1,
                         labelWidth: options.labelWidth,
                         ticks: [0, maxInt*0.1, maxInt*0.2, maxInt*0.3, maxInt*0.4, maxInt*0.5,
                                 maxInt*0.6, maxInt*0.7, maxInt*0.8, maxInt*0.9, maxInt],
                         tickFormatter: function(val, axis) {return Math.round((val * 100)/maxInt)+"%";}}
            };
        container.data("plotOptions", plotOptions);
        container.data("maxInt", maxInt);

    }

    function getMaxInt(options) {
        var maxInt = 0;
        for(var j = 0; j < options.spectrum.length; j += 1) {
            var peak = options.spectrum[j];
            if(peak[1] > maxInt) {
                maxInt = peak[1];
            }
        }
        if (options.spectrum.length == 1)
            maxInt = maxInt * 1.5;
        //alert(maxInt);
        return maxInt;
    }

    function setMassError(container) {
        var me = parseFloat($(getElementSelector(container, elementIds.massError)).val());
        var unit = getMassErrorUnit(container);
        if(me != container.data("massError")) {
            container.data("massError", me);
            container.data("massErrorChanged", true);
        }
        else if(options.massErrorUnit !== unit)
        {
            options.massErrorUnit = unit;
            container.data("massErrorChanged", true);
        }
		else {
			container.data("massErrorChanged", false);
		}
    }

    // -----------------------------------------------
    // CREATE MS/MS PLOT
    // -----------------------------------------------
    function createPlot(container, datasets) {

        var plot;
        var options = container.data("options");
        if(!container.data("zoomRange")) {
        // Set the default X range to be from 50 to the MW of the peptide.
        // This allows easier comparison of different spectra from the
                // same peptide ion in different browser tabs. One can blink back
                // and forth when the display range is identical.
            var selectOpts = {};
            //var neutralMass = options.peptide.getNeutralMassMono() + Ion.MASS_O + Ion.MASS_H;
            plot = $.plot(getElementSelector(container, elementIds.msmsplot), datasets,
                      $.extend(true, {}, container.data("plotOptions"), selectOpts));
            //plot = $.plot($(getElementSelector(container, elementIds.msmsplot)), datasets,  container.data("plotOptions"));
        }
        else {
            var zoomRange = container.data("zoomRange");
            var selectOpts = {};
            if($(getElementSelector(container, elementIds.zoom_x)).is(":checked"))
                selectOpts.xaxis = { min: zoomRange.xaxis.from, max: zoomRange.xaxis.to };
            if($(getElementSelector(container, elementIds.zoom_y)).is(":checked"))
                selectOpts.yaxis = { min: 0, max: zoomRange.yaxis.to };

            plot = $.plot(getElementSelector(container, elementIds.msmsplot), datasets,
                      $.extend(true, {}, container.data("plotOptions"), selectOpts));

            // zoom out icon on plot right hand corner
            var o = plot.getPlotOffset();
            $(getElementSelector(container, elementIds.msmsplot)).append('<div id="'+getElementId(container, elementIds.ms2plot_zoom_out)+'" class="zoom_out_link" style="position:absolute; left:'
                    + (o.left + plot.width() - 20) + 'px;top:' + (o.top+4) + 'px"></div>');

            $(getElementSelector(container, elementIds.ms2plot_zoom_out)).click( function() {
                resetZoom(container, options);
            });
        }

        // we have re-calculated and re-drawn everything..
        container.data("massTypeChanged", false);
        container.data("massErrorChanged", false);
        container.data("plot", plot);
        // Draw the peak mass error plots
        plotPeakMassErrorPlot(container, datasets);
        if(container.data("options").showMassErrorPlot === false)
        {
            $(getElementSelector(container, elementIds.massErrorPlot)).hide();
        }
    }

    function getPlotXRange(options) {
        var xmin = options.spectrum[0][0];
        var xmax = options.spectrum[options.spectrum.length - 1][0];
        var xpadding = (xmax - xmin) * 0.025;
        if (xpadding == 0) xpadding = 10;
        return {xmin:xmin - xpadding, xmax:xmax + xpadding};
    }

    function plotPeakMassErrorPlot(container, datasets) {
        var data = [];
        var options = container.data("options");

        var minMassError = 0;
        var maxMassError = 0;

        var s_data = [];
        for (var i = 0; i < options.annotation.length; i += 1) {
            var massError = options.massErrors[i];
            var mz = options.matchedPeaks[i][0];
            minMassError = Math.min(minMassError, massError);
            maxMassError = Math.max(maxMassError, massError);
            s_data.push([mz, massError]);
        }

        var placeholder = $(getElementSelector(container, elementIds.massErrorPlot));

        var __xrange = getPlotXRange(options);
        var zoomRange = container.data("zoomRange");
        if (zoomRange)
        {
            // Sync zooming with the MS/MS plot.
            __xrange.xmin = zoomRange.xaxis.from;
            __xrange.xmax = zoomRange.xaxis.to;
        }

        var balancedMassError = Math.max(Math.abs(maxMassError), Math.abs(minMassError)) * 1.05;
        var tickY = (balancedMassError / 2).toFixed(2);

        // the MS/MS plot should have been created by now.  This is a hack to get the plots aligned.
        // We will set the y-axis labelWidth to this value.
        var labelWidth = container.data("plot").getAxes().yaxis.labelWidth;

        var massErrorPlotOptions = {
            series: {data: s_data},
            grid: {
                show: true,
                autoHighlight: false,
                clickable: false,
                hoverable: true,
                mouseActiveRadius: 5,
                borderWidth: 1,
                labelMargin: options.labelMargin,
                markings: [
                    {yaxis: {from: 0, to: 0}, color: "#555555", lineWidth: 0.5}
                ]  // draw a horizontal line at y=0
            },
            selection: {mode: "xy", color: "#F0E68C"},
            xaxis: {
                tickLength: 3, tickColor: "#000",
                min: __xrange.xmin,
                max: __xrange.xmax
            },
            yaxis: {
                tickLength: 3, tickColor: "#000",
                ticks: [-tickY, 0, tickY],
                min: -balancedMassError,
                max: balancedMassError,
                labelWidth: labelWidth
            }
        };
        data.push({data: s_data,
                  color: "#000",
                  points: {show: true, fill: true, radius: 1},
                  labelType: 'none'});

        // TOOLTIPS
        $(getElementSelector(container, elementIds.massErrorPlot)).bind("plothover", function (event, pos, item) {
            displayTooltip(item, container, options, "M/z", "Error");
        });

        // CHECKBOX TO SHOW OR HIDE MASS ERROR PLOT
        var massErrorPlot = $.plot(placeholder, data, massErrorPlotOptions);

        // Display clickable mass error unit.
        // TODO: if we need this, we should make it work with jquery-1.4.2 (see Lorikeet for reference)
//        var o = massErrorPlot.getPlotOffset();
//        placeholder.append('<div id="' + getElementId(container, elementIds.massErrorPlot_unit) + '" class="link"  '
//            + 'style="position:absolute; left:'
//            + (o.left + 5) + 'px;top:' + (o.top + 4) + 'px;'
//            + 'background-color:yellow; font-style:italic">'
//            + options.massErrorUnit + '</div>');
//        $(getElementSelector(container, elementIds.massErrorPlot_unit)).click(function () {
//            var unit = $(this).text();
//
//            if (unit == massErrorType_Da) {
//                $(this).text(massErrorType_ppm);
//                options.massErrorUnit = massErrorType_ppm;
//                $("input[name='"+getRadioName(container, "massErrorUnit")+"'][value='" + massErrorType_ppm + "']").prop("checked",true)
//            }
//            else if (unit == massErrorType_ppm) {
//                $(this).text(massErrorType_Da);
//                options.massErrorUnit = massErrorType_Da;
//                $("input[name='"+getRadioName(container, "massErrorUnit")+"'][value='" + massErrorType_Da + "']").prop("checked",true)
//            }
//
//            changeMassErrorUnit(options, container, datasets);
//        });
        var o = massErrorPlot.getPlotOffset();
    }

    function changeMassErrorUnit(options, container, datasets) {
        refreshMassErrors(options);
        plotPeakMassErrorPlot(container, datasets);
        addMatchedPeaksTable(container);
        if (options.selectedPeak) {
            selectGroup(container, options.selectedPeak, options);
        }
    }

    function displayTooltip(item, container, options, tooltip_xlabel, tooltip_ylabel) {
        if (item) {
            if (container.data("previousPoint") != item.datapoint) {
                container.data("previousPoint", item.datapoint);

                $(getElementSelector(container, elementIds.msmstooltip)).remove();
                var x = item.datapoint[0].toFixed(options.precision),
                    y = item.datapoint[1].toFixed(getMassErrorPrecision(options));

                showTooltip(container, item.pageX, item.pageY,
                    tooltip_xlabel + ": " + x + "<br>" + tooltip_ylabel + ": " + y, options);
            }
        }
        else {
            $(getElementSelector(container, elementIds.msmstooltip)).remove();
            container.data("previousPoint", null);
        }
    }

    // -----------------------------------------------
    // SET UP INTERACTIVE ACTIONS FOR MS/MS PLOT
    // -----------------------------------------------
    function setupInteractions (container, options) {

        // ZOOMING
        $(getElementSelector(container, elementIds.msmsplot)).bind("plotselected", function (event, ranges) {
            container.data("zoomRange", ranges);
            reloadPlotWithSelection(container, options);
        });

        // ZOOM AXES
        $(getElementSelector(container, elementIds.zoom_x)).click(function() {
            resetAxisZoom(container);
        });
        $(getElementSelector(container, elementIds.zoom_y)).click(function() {
            resetAxisZoom(container);
        });

        // RESET ZOOM
        $(getElementSelector(container, elementIds.resetZoom)).click(function() {
            resetZoom(container, options);
        });

        // UPDATE
        $(getElementSelector(container, elementIds.update)).click(function() {
            container.data("zoomRange", null); // zoom out fully
            setMassError(container);
            createPlot(container, getDatasets(container));
        });

        $(getElementSelector(container, elementIds.msmsplot)).bind("plotclick", function (event, pos, item) {
                var plotOptions = container.data("options");
                if (item) {
                    var selectedPeakIdx = -1;
                    for (var i = 0; i < plotOptions.matchedPeaks.length; i++) {
                        if(item.datapoint[0] == plotOptions.matchedPeaks[i][0]) {
                            selectedPeakIdx = i;
                            break;
                        }
                    }
                    if (selectedPeakIdx != -1){
                        selectGroup(container, selectedPeakIdx, plotOptions);
                        container.data("options", plotOptions);
                        createPlot(container, getDatasets(container, plotOptions.matchedPeaks[selectedPeakIdx]));
                    }
                }
            });
        $(getElementSelector(container, elementIds.enableTooltip)).click(function() {
            $(getElementSelector(container, elementIds.msmstooltip)).remove();
        });

        // PLOT MASS ERROR CHECKBOX
        $(getElementSelector(container, elementIds.massErrorPlot_option)).click(function() {
            var plotDiv = $(getElementSelector(container, elementIds.massErrorPlot));
            if($(this).is(':checked'))
            {
                plotDiv.show();
                plotPeakMassErrorPlot(container, getDatasets(container));
                container.data("options").showMassErrorPlot = true;
            }
            else
            {
                plotDiv.hide();
                container.data("options").showMassErrorPlot = false;
            }
        });

        // CHANGING THE PLOT SIZE
        makePlotResizable(container);

	    // PRINT SPECTRUM
	    savePlot(container);
    }

    function isComponentModified(options, atomComponentId) {
        if (!options.modification)
            return false;
        var modifiedComponent = options.modification.position;
        return atomComponentId == modifiedComponent;
    }

    function showNewMol(mol, options) {
        var fragmentedBonds = options.fragmentedBonds;
        for(var i = 0; i < mol.atoms.length; i++) {
            var atomComponentId = options.atomComponents[i];
            if (isComponentModified(options, atomComponentId)) {
                mol.atoms[i].specs = modifiedSpecs;
            }
            else mol.atoms[i].specs = nonSelectedSpecs;
        }
        for(var i = 0; i < mol.bonds.length; i++) {
            var b = mol.bonds[i];
            if (fragmentedBonds.indexOf(i) != -1) {
                b.specs = fragmentedBondsSpecs;
            }
            else if (b.a1.specs == modifiedSpecs && b.a2.specs == modifiedSpecs) {
                b.specs = modifiedSpecs;
            }
            else b.specs = nonSelectedSpecs;
        }
        options.canvas.loadMolecule(mol);
        return mol;
    }

    function resetSelection(container) {
        var options = container.data("options");

        if (options.selectedPeak != null) {
            var matchedPeadId = getElementId(container, elementIds.peakRow) + '_' + options.selectedPeak;
            $('#' + matchedPeadId).removeClass('selected_row');
            options.selectedPeak = null;

            if (options.mol) {
                showNewMol(options.mol, options);
            }
            else {
                $(getElementSelector(container, elementIds.seqCanvas)).html(getModifiedSequence(options));
            }
            createPlot(container, getDatasets(container));
            $(getElementSelector(container, elementIds.currentMass))
                .text('Select a colored peak to view annotation');
        }
    }

    function selectGroup(container, matchedPeakIdx, options) {
        options.selectedPeak = matchedPeakIdx;
        var mz = options.matchedPeaks[matchedPeakIdx][0];
        var intensity = options.matchedPeaks[matchedPeakIdx][1];
        var massError = options.massErrors[matchedPeakIdx];
        var charge = options.annotation[matchedPeakIdx].charge;
        var components = options.annotation[matchedPeakIdx].components;

        // clearing selection for other peaks
        for (var i = 0; i < options.matchedPeaks.length; i++) {
            $('#' + getElementId(container, elementIds.peakRow) + '_' + i).attr('class', '');
        }

        var matchInfo = [];
        matchInfo.push(parseFloat(mz).toFixed(options.precision));
        matchInfo.push(parseFloat(massError).toFixed(getMassErrorPrecision(options)));
        matchInfo.push(charge);
        matchInfo.push(parseFloat(intensity).toFixed(options.precision));
        var matchedPeadId = getElementId(container, elementIds.peakRow) + '_' + matchedPeakIdx;
        $('#' + matchedPeadId).addClass('selected_row');
        //document.getElementById(matchedPeadId).scrollIntoView(true);

        if (options.mol)
            selectGroupInMol(components, options);
        else
            selectGroupInSeq(container, components, options);
        $(getElementSelector(container, elementIds.currentMass)).text(
            matchInfoToString(matchInfo, options.massErrorUnit));
    }

    function matchInfoToString(matchInfo, unit) {
        return "M/z: " + matchInfo[0] + " Da/e, " + "mass error: " + matchInfo[1] +
            " " + unit + ", charge: " + matchInfo[2] + ", intensity: " + matchInfo[3];
    }

    function selectGroupInSeq(container, components, options) {
        $(getElementSelector(container, elementIds.seqCanvas)).html(getModifiedSequence(options, components));
    }

    function selectGroupInMol(components, options) {
        var mol = options.mol;
        var fragmentedBonds = options.fragmentedBonds;
        for (var i = 0; i < mol.atoms.length; i++) {
            var atomComponentId = options.atomComponents[i];
            var isSelected = components.indexOf(atomComponentId) != -1;
            var isModified = isComponentModified(options, atomComponentId);
            if (isModified && isSelected) {
                mol.atoms[i].specs = modSelectedSpecs;
            }
            else if (isModified) {
                mol.atoms[i].specs = modifiedSpecs;
            }
            else if (isSelected) {
                mol.atoms[i].specs = selectedSpecs;
            }
            else mol.atoms[i].specs = nonSelectedSpecs;
        }
        for (var i = 0; i < mol.bonds.length; i++) {
            var b = mol.bonds[i];
            if (b.a1.specs == invisibleSpecs || b.a2.specs == invisibleSpecs)
                b.specs = invisibleSpecs;
            else if (fragmentedBonds.indexOf(i) != -1)
                b.specs = fragmentedBondsSpecs;
            else if (b.a1.specs == modifiedSpecs && b.a2.specs == modifiedSpecs)
                b.specs = modifiedSpecs;
            else if (b.a1.specs == modSelectedSpecs && b.a2.specs == modSelectedSpecs)
                b.specs = modSelectedSpecs;
            else if (b.a1.specs == selectedSpecs && b.a2.specs == selectedSpecs)
                b.specs = selectedSpecs;
            else b.specs = nonSelectedSpecs;
        }
        options.canvas.loadMolecule(mol);
    }

    function addMatchedPeaksTable(container) {
        var options = container.data("options");

        var matchInfos = [];
        for (var matchedPeakIdx = 0; matchedPeakIdx < options.annotation.length; matchedPeakIdx++) {
            var mz = options.matchedPeaks[matchedPeakIdx][0];
            var intensity = options.matchedPeaks[matchedPeakIdx][1];
            var massError = options.massErrors[matchedPeakIdx];
            var charge = options.annotation[matchedPeakIdx].charge;

            var matchInfo = [];
            matchInfo.push(parseFloat(mz).toFixed(options.precision));
            matchInfo.push(parseFloat(massError).toFixed(getMassErrorPrecision(options)));
            matchInfo.push(charge);
            matchInfo.push(parseFloat(intensity).toFixed(options.precision));
            matchInfos.push(matchInfo);
        }

        var peaksTable = '' ;
        peaksTable += '<div id="' + getElementId(container, elementIds.peaksTable) + '" align="center">';
        peaksTable += '<h4 class="annoPeaksHeader">Annotated peaks</h4>';
        peaksTable += '<table id="table_scroll" cellpadding="1" class="font_small ' + 'annoPeaksTable' + '">';
        peaksTable +=  '<thead>';
        peaksTable +=   "<tr>";
        var headersPeaksTable = ["M/z (Da/e)", "Mass error (" + options.massErrorUnit + ")",
            "Charge", "Intensity"];

        for(var i = 0; i < headersPeaksTable.length; i += 1) {
            peaksTable += '<th>' + headersPeaksTable[i] +  '</th>';
        }
        peaksTable += "</tr>";
        peaksTable += "</thead>";

        if (options.mol)
            peaksTable += '<tbody style="height: 170px;">';
        else
            peaksTable += '<tbody style="height: 370px;">';

        for(var i = 0; i < matchInfos.length; i += 1) {
            peaksTable +=   '<tr id="' + getElementId(container, elementIds.peakRow) + '_' + i + '">';
            for (var j = 0; j < matchInfos[i].length; j++){
                peaksTable += "<td>" + matchInfos[i][j] +  "</td>";
            }
            peaksTable += "</tr>";
        }

        peaksTable += "</tbody>";
        peaksTable += "</table>";
        peaksTable += "</div>";

        $(getElementSelector(container, elementIds.peaksTable)).remove();
        $(getElementSelector(container, elementIds.peaksTableDiv)).prepend(peaksTable);

        if ( options.sizeChangeCallbackFunction ) {
            options.sizeChangeCallbackFunction();
        }
    }

    function resetZoom(container, options) {
        container.data("zoomRange", null);
        setMassError(container);
        reloadPlotWithSelection(container, options);
    }

    function reloadPlotWithSelection(container, options) {
        if (options.selectedPeak) {
            plotOptions = container.data("options");
            createPlot(container, getDatasets(container, plotOptions.matchedPeaks[options.selectedPeak]));
        }
        else {
            createPlot(container, getDatasets(container));
        }
    }

    function plotAccordingToChoices(container) {
        var data = getDatasets(container);

        if (data.length > 0) {
            createPlot(container, data);
            showSequenceInfo(container); // update the MH+ and m/z values
        }
    }

    function resetAxisZoom(container) {

        var plot = container.data("plot");
        var plotOptions = container.data("plotOptions");

        var zoom_x = false;
        var zoom_y = false;
        if($(getElementSelector(container, elementIds.zoom_x)).is(":checked"))
            zoom_x = true;
        if($(getElementSelector(container, elementIds.zoom_y)).is(":checked"))
            zoom_y = true;
        if(zoom_x && zoom_y) {
            plotOptions.selection.mode = "xy";
            if(plot) plot.getOptions().selection.mode = "xy";
        }
        else if(zoom_x) {
            plotOptions.selection.mode = "x";
            if(plot) plot.getOptions().selection.mode = "x";
        }
        else if(zoom_y) {
            plotOptions.selection.mode = "y";
            if(plot) plot.getOptions().selection.mode = "y";
        }
    }

    function showTooltip(container, x, y, contents, options) {

        var tooltipCSS = {
            position: 'absolute',
            display: 'none',
            top: y + 5,
            left: x + 5,
            border: '1px solid #fdd',
            padding: '2px',
            'background-color': '#F0E68C',
            opacity: 1 };

        $('<div id="'+getElementId(container, elementIds.msmstooltip)+'">' + contents + '</div>')
                .css( tooltipCSS ).appendTo("body").fadeIn(200);
    }

    function makePlotResizable(container) {

        var options = container.data("options");

        $(getElementSelector(container, elementIds.slider_width)).slider({
            value:options.width,
            min: 100,
            max: 1500,
            step: 50,
            slide: function(event, ui) {
                var width = ui.value;
                options.width = width;
                $(getElementSelector(container, elementIds.msmsplot)).css({width: width});
                $(getElementSelector(container, elementIds.massErrorPlot)).css({width: width});

                plotAccordingToChoices(container);
                $(getElementSelector(container, elementIds.slider_width_val)).text(width);
                if ( options.sizeChangeCallbackFunction ) {
                    options.sizeChangeCallbackFunction();
                }
            }
        });

        $(getElementSelector(container, elementIds.slider_height)).slider({
            value:options.height,
            min: 100,
            max: 1000,
            step: 50,
            slide: function(event, ui) {
                var height = ui.value;
                options.height = height;
                $(getElementSelector(container, elementIds.msmsplot)).css({height: height});
                plotAccordingToChoices(container);
                $(getElementSelector(container, elementIds.slider_height_val)).text(height);
                if ( options.sizeChangeCallbackFunction ) {
                    options.sizeChangeCallbackFunction();
                }
            }
        });
    }

	function savePlot(container) {
		$(getElementSelector(container, elementIds.saveLink)).click(function() {
            var parent = container.parent();

			// create another div and move the plots into that div
			$(document.body).append('<div id="tempPrintDiv"></div>');
			$("#tempPrintDiv").append(container.detach());
			$("#tempPrintDiv").siblings().addClass("noprint");

			var plotOptions = container.data("plotOptions");

			container.find(".bar").addClass('noprint');
			$(getElementSelector(container, elementIds.peaksTableDiv)).addClass('noprint');
			$(getElementSelector(container, elementIds.chemCanvas)).parent().addClass('noprint');

            var plotCanvas = container.data("plot").getCanvas();
            var plotContext = plotCanvas.getContext("2d");
            plotContext.globalCompositeOperation = "destination-over";
            plotContext.fillStyle = '#fff';
            plotContext.fillRect(0, 0, plotCanvas.width, plotCanvas.height);

			var image = plotCanvas.toDataURL("image/png").replace("image/png", "image/octet-stream");
            var link = document.createElement('a');
            link.href = image;
            link.download = location.pathname.substring(location.pathname.lastIndexOf("/") + 1).replace('.html', '.png');
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
			// remove the class after printing so that if the user prints
			// via the browser's print menu the whole page is printed
			container.find(".bar").removeClass('noprint');
			$(getElementSelector(container, elementIds.peaksTableDiv)).removeClass('noprint');
			$(getElementSelector(container, elementIds.chemCanvas)).parent().removeClass('noprint');
			$("#tempPrintDiv").siblings().removeClass("noprint");

			// move the plots back to the original location
            parent.append(container.detach());
			$("#tempPrintDiv").remove();

		});
	}

    function getMassErrorUnit(container)
    {
        return container.find("input[name='"+getRadioName(container, "massErrorUnit")+"']:checked").val();
    }

    // -----------------------------------------------
    // SELECTED DATASETS
    // -----------------------------------------------
    function getDatasets(container, selectedPeak) {

        var options = container.data("options");

        var data = [{data: options.notMatchedPeaks, color: '#C0C0C0', bars: {
                        show: true,
                        fill: 1,
                        barWidth: 2,
                        lineWidth: 0.5,
                        fillColor:  "#C0C0C0",
                        align: 'center'
                    }, clickable: false, hoverable: false},
                    {data: options.matchedPeaks, color: '#00CCFF', bars: {
                        show: true,
                        fill: 1,
                        lineWidth: 0.5,
                        color: '#66CCFF',
                        barWidth: 3,
                        align: 'center'
                    }}];
        if (selectedPeak){
            var d = {data: [selectedPeak], color: '#0000FF', bars: {
                        show: true,
                        fill: 1,
                        lineWidth: 0.5,
                        color: '#0000FF',
                        barWidth: 3,
                        align: 'center'}};
            data.push(d);
        }

        return data;
    }

    function processInfoField(field, unit, precision) {
        if (unit)
            unit = ' ' + unit;
        else
            unit = '';
        return (field || field === 0) ? (precision ? field.toFixed(precision) : field) + unit : 'N/A';
    }

    // -----------------------------------------------
    // INITIALIZE THE CONTAINER
    // -----------------------------------------------
    function addGeneralInfoPanel(div, options) {
        div.append('<div id="'+elementIds.psm_general_info+'"></div>');
        var infoContainer = $("#"+ div.attr('id')+" > #"+elementIds.psm_general_info);

        if (options.spectrumInfo) {
            var table = '<table class="infoTable"> ';
            table += '<thead>';
            table += '<tr>';
            table += '<th colspan="2">Spectrum</th>';
            table += '</tr></thead>';
            table += '<tbody> ';
            table += '<tr>';
            table += '<td class="field_name"><span>Filename</span></td>';
            table += '<td class="field_value"><span>' + trimName(processInfoField(options.spectrumInfo.filename)) + '</span></td>';
            table += '</tr><tr>';
            table += '<td class="field_name"><span>Scan</span></td>';
            table += '<td class="field_value"><span>' + processInfoField(options.spectrumInfo.scan) + '</span></td>';
            table += '</tr><tr>';
            table += '<td class="field_name"><span>M/z</span></td>';
            table += '<td class="field_value"><span>' + processInfoField(options.spectrumInfo.mz, 'Da/e', options.precision) + '</span></td>';
            table += '</tr><tr>';
            table += '<td class="field_name"><span>Charge</span></td>';
            table += '<td class="field_value"><span>' + processInfoField(options.spectrumInfo.charge) + '</span></td>';
            table += '</tr><tr>';
            table += '<td class="field_name"><span>Retention</span></td>';
            table += '<td class="field_value"><span>' + processInfoField(options.spectrumInfo.retention, 'sec', options.precision) + '</span></td>';
            table += '<tr>';
            table += '</table>';
            infoContainer.append(table);
        }

        if (options.compoundInfo) {
            var table = '<table class="infoTable"> ';
            table += '<thead>';
            table += '<tr>';
            table += '<th colspan="2">Compound</th>';
            table += '</tr></thead>';
            table += '<tbody> ';
            table += '<tr>';
            table += '<td class="field_name"><span>Name</span></td>';
            table += '<td class="field_value"><span>' + trimName(processInfoField(options.compoundInfo.name)) + '</span></td>';
            table += '</tr><tr>';
            table += '<td class="field_name"><span>Mass</span></td>';
            table += '<td class="field_value"><span>' + processInfoField(options.compoundInfo.mass, 'Da', options.precision) + '</span></td>';
            table += '</tr><tr>';
            table += '<td class="field_name"><span>Mass error</span></td>';
            table += '<td class="field_value"><span>' + processInfoField(options.compoundInfo.massError, 'Da', options.precision) + '</span></td>';
            table += '</tr><tr>';
            table += '<td class="field_name"><span>Formula</span></td>';
            table += '<td class="field_value"><span>' + processInfoField(options.compoundInfo.formula) + '</span></td>';
            table += '</tr><tr>';
            table += '<td class="field_name"><span>Adduct</span></td>';
            table += '<td class="field_value"><span>' + processInfoField(options.compoundInfo.adduct) + '</span></td>';
            table += '</tr>';
            table += '</table>';
            infoContainer.append(table);
        }
    }

    function trimName(name) {
        var maxLength = 50;
        if (name.length <= maxLength)
            return name;

        var separator = '...';
        var charsToShow = maxLength - separator.length,
            firstChars = Math.ceil(charsToShow / 2),
            lastChars = Math.floor(charsToShow / 2);

        return name.substr(0, firstChars) +
               separator +
               name.substr(name.length - lastChars);
    }

    function createContainer(div) {
        var contentId = elementIds.lorikeet_content + "_" + index;
        div.append('<div id="'+ contentId +'" tabindex="' + index + '"></div>');
        var container = $("#"+ div.attr('id')+" > #" + contentId);
        container.addClass("lorikeet");
        $("#" + contentId).keyup(function(e) {
            if (e.keyCode === 27) resetSelection(container);   // esc
        });
        return container;
    }

    function initContainer(container) {

        var options = container.data("options");

        var rowspan = 2;

        var parentTable = '<table cellpadding="0" cellspacing="5" class="lorikeet-outer-table"> ';
        parentTable += '<tbody> ';

        if(options.sequence) {
            // placeholder for sequence, m/z, scan number etc
            parentTable += '<td colspan="2" style="background-color: white; padding:5px; border:1px dotted #cccccc;" valign="bottom" align="center"> ';
            parentTable += '<div id="'+getElementId(container, elementIds.seqinfo)+'" style="width:100%;"></div> ';
            parentTable += '</td> ';
        }

        // placeholders for the ms/ms plot
        parentTable += '<tr> ';
        parentTable += '<td style="background-color: white; padding:5px; border:1px dotted #cccccc;width:'+options.width+'px;height:'+options.height+'px;" valign="top" align="center"> ';
        // placeholder for peak mass error plot
        parentTable += '<div id="'+getElementId(container, elementIds.msmsplot)+'" align="bottom" style="width:'+options.width+'px;height:'+options.height+'px;"></div> ';
        parentTable += '<div id="'+getElementId(container, elementIds.viewOptionsDiv)+'" style="margin-top:15px;"></div> ';
        parentTable += '<div id="'+getElementId(container, elementIds.currentMass)+'" style="margin-top:15px;">Select a colored peak to view annotation</div> ';
        // placeholder for viewing options (zoom, plot size etc.)
        parentTable += '<div id="'+getElementId(container, elementIds.massErrorPlot)+'" style="width:'+options.width+'px;height:90px;"></div> ';
        parentTable += '</td> ';
        parentTable += '<td style="background-color: white; padding:5px; border:1px dotted #cccccc;width:" valign="top" align="center"> ';
        if (options.modification) {
            parentTable += modificationInfoToString(options.modification, options.precision);
        }
        parentTable += '<canvas id="' + getElementId(container, elementIds.chemCanvas) + '" align="top"></canvas>';
        if (options.mol) {
            parentTable += '<div style="margin-bottom: 20px;"><span class="span_tip">Click and use mouse to roll and zoom the structure</span></div>';
        }

        parentTable += '<div id="'+ getElementId(container, elementIds.peaksTableDiv) +'" align="top"></div> ';
        parentTable += '</td> ';
        parentTable += '</tr> ';

        parentTable += '</tbody> ';
        parentTable += '</table> ';

        container.append(parentTable);
        if (options.mol) {
            var myCanvas = new ChemDoodle.TransformCanvas('' + getElementId(container, elementIds.chemCanvas), 550, 250);
            myCanvas.emptyMessage = 'No Data Loaded!';
            myCanvas.loadMolecule(options.mol);
            options.canvas = myCanvas;
        }
        else {
            $('#' + getElementId(container, elementIds.chemCanvas)).hide();
        }
        return container;
    }

    //---------------------------------------------------------
    // SEQUENCE INFO
    //---------------------------------------------------------
    function showSequenceInfo (container) {
        var options = container.data("options");

        var specinfo = '';
        if(options.sequence) {
            specinfo += '<div>';
            specinfo += '<p style="font-weight:bold; color:#8B0000;" id="' + getElementId(container, elementIds.seqCanvas) + '">' + getModifiedSequence(options) + '</p>';
            specinfo += '</div>';
        }

        // first clear the div if it has anything
        $(getElementSelector(container, elementIds.seqinfo)).empty();
        $(getElementSelector(container, elementIds.seqinfo)).append(specinfo);
    }

    function getModifiedSequence(options, proteinArray) {
        var modSeq = '';
        var lettersCounter = 0;
        var altLetters = [];
        var altPattern = /\w[-+]\d/ig;
        while (result = altPattern.exec(options.sequence)) {
            altLetters.push(result.index);
        }
        for(var i = 0; i < options.sequence.length; i += 1) {
            currentChar = options.sequence.charAt(i);
            if (currentChar.match(/[A-Z]/i)) lettersCounter += 1;
            isAltLetter = altLetters.indexOf(i) != -1 || !currentChar.match(/[A-Z]/i);
            isSelected = proteinArray && proteinArray.indexOf(lettersCounter - 1) != -1;
            if (isAltLetter && isSelected)
                modSeq += highlightSequence(currentChar, '#00ccff');
            else if (isSelected)
                modSeq += highlightSequence(currentChar, 'blue');
            else if (isAltLetter)
                modSeq += highlightSequence(currentChar, 'red');
            else modSeq += currentChar;
        }
        return modSeq;
    }

    function highlightSequence(seq, color) {
        return '<span style="color: ' + color + ';">' + seq + '</span>';
    }

    function modificationInfoToString(modification, precision) {
        var massShiftSign = modification.massShift > 0 ? ' + ' : ' - ';
        var massShiftStr = massShiftSign + Math.abs(modification.massShift).toFixed(precision);
        var componentActualMass = Math.max(0, modification.componentMass + modification.massShift).toFixed(precision);
        var lossAA = componentActualMass < Ion.MASS_PROTON ? ' (loss of AA)' : '';
        return '<div id="mass_shift_div" align="left" style="padding-left: 20px; padding-top: 10px;z-index: 1000;">' +
                '<span style="font-weight:bold;">' +
                'Modification: <span style="color: ' + molColors['modified'] + ';">' +
                modification.componentFormula + massShiftStr + ' Da = ' +
                componentActualMass + ' Da</span>' + lossAA + '</div>';
    }

    //---------------------------------------------------------
    // VIEWING OPTIONS TABLE
    //---------------------------------------------------------
    function makeViewingOptions(container) {
        var options = container.data("options");

        var myContent = '';

        // reset zoom option
        myContent += '<nobr> ';
        myContent += '<span class="span_tip">Click and drag in the plot to zoom</span> ';
        myContent += 'X:<input id="'+getElementId(container, elementIds.zoom_x)+'" type="checkbox" value="X" checked="checked"/> ';
        myContent += '&nbsp;Y:<input id="'+getElementId(container, elementIds.zoom_y)+'" type="checkbox" value="Y" /> ';
        myContent += '&nbsp;<input id="'+getElementId(container, elementIds.resetZoom)+'" type="button" value="Zoom Out" /> ';
		myContent += '&nbsp;<input id="'+getElementId(container, elementIds.saveLink)+'" type="button" value="Save as PNG" /> ';
        myContent += '</nobr> ';

        myContent += '&nbsp;&nbsp;';

        // mass error plot option
        myContent += '<br>';
        myContent += 'Mass errors: ' +
            '<label><input type="radio" name="'+getRadioName(container, "massErrorUnit")+'" value="' + massErrorType_Da + '"';
        if(options.massErrorUnit === massErrorType_Da)
        {
            myContent += ' checked = "checked" ';
        }
        myContent += '/><span style="font-weight: bold;">' + massErrorType_Da + '</span></label> ';
        myContent += '<label><input type="radio" name="'+getRadioName(container, "massErrorUnit")+'" value="' + massErrorType_ppm + '"';
        if(options.massErrorUnit === massErrorType_ppm)
        {
            myContent += ' checked = "checked" ';
        }
        myContent += '/><span style="font-weight: bold;">' + massErrorType_ppm + '</span></label> ';
        myContent += '<label><input id="'+getElementId(container, elementIds.massErrorPlot_option)+'" type="checkbox" ';
        if(options.showMassErrorPlot === true)
        {
            myContent += 'checked="checked"';
        }
        myContent += '> Plot</label> ';
        $(getElementSelector(container, elementIds.viewOptionsDiv)).append(myContent);
        if(!options.showViewingOptions) {
            $(getElementSelector(container, elementIds.viewOptionsDiv)).hide();
        }

        $("input[name='"+getRadioName(container, "massErrorUnit")+"']").click(function () {
            var unit = $(this).val();

            if (unit == massErrorType_ppm) {
                $(getElementSelector(container, elementIds.massErrorPlot_unit)).text(massErrorType_ppm);
                options.massErrorUnit = massErrorType_ppm;
            }
            else if (unit == massErrorType_Da) {
                $(getElementSelector(container, elementIds.massErrorPlot_unit)).text(massErrorType_Da);
                options.massErrorUnit = massErrorType_Da;
            }

            changeMassErrorUnit(options, container, getDatasets(container));
        });

    }

});
