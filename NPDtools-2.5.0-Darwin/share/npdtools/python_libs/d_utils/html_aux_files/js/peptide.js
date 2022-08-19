// $LastChangedDate$
// $LastChangedBy$
// $LastChangedRevision$

// -----------------------------------------------------------------------------
// Peptide sequence and modifications
// -----------------------------------------------------------------------------
function Peptide(seq, staticModifications, varModifications, ntermModification, ctermModification, maxNeutralLossCount) {
	
	var sequence = seq;
    if(!sequence) {
        sequence = "";
    }

	var ntermMod = ntermModification;
	var ctermMod = ctermModification;
	var staticMods = [];
    var varMods = [];
    var potentialLosses_custom = {};
    var potentialLosses_lorikeet = {};
    var potentialLossesAtIndex = [];
    var nterm_totalLossOptions = [];
    var cterm_totalLossOptions = [];
    var maxNeutralLossCount = maxNeutralLossCount;

    var debug = false;

    _initNeutralLosses();
    _initStaticMods();
    _initVarMods();
    // calculate the total loss options at each index using only custom neutral losses.
    _calculateTotalLossOptions(null, maxNeutralLossCount);


    //----------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------
    this.sequence = function () {
        return sequence;
    }

    this.varMods = function() {
        return varMods;
    }

    // index: index in the seq.
    // If this is a N-term sequence we will sum up the mass of the amino acids in the sequence up-to index (exclusive).
    // If this is a C-term sequence we will sum up the mass of the amino acids in the sequence starting from index (inclusive)
    // modification masses are added
    this.getSeqMassMono = function _seqMassMono(index, term) {
        return _getSeqMass(index, term, "mono");
    }

    // index: index in the seq.
    // If this is a N-term sequence we will sum up the mass of the amino acids in the sequence up-to index (exclusive).
    // If this is a C-term sequence we will sum up the mass of the amino acids in the sequence starting from index (inclusive)
    // modification masses are added
    this.getSeqMassAvg = function _seqMassAvg(index, term) {
        return _getSeqMass(index, term, "avg");
    }

    // Returns the monoisotopic neutral mass of the peptide; modifications added. N-term H and C-term OH are added
    this.getNeutralMassMono = function _massNeutralMono() {

        var mass = 0;
        var aa_obj = new AminoAcid();
        if(sequence) {
            for(var i = 0; i < sequence.length; i++) {
                var aa = aa_obj.get(sequence.charAt(i));
                mass += aa.mono;
            }
        }

        mass = _addTerminalModMass(mass, "n");
        mass = _addTerminalModMass(mass, "c");
        mass = _addResidueModMasses(mass, sequence.length, "n");
        // add N-terminal H
        mass = mass + Ion.MASS_H_1;
        // add C-terminal OH
        mass = mass + Ion.MASS_O_16 + Ion.MASS_H_1;

        return mass;
    }

    //Returns the avg neutral mass of the peptide; modifications added. N-term H and C-term OH are added
    this.getNeutralMassAvg = function _massNeutralAvg() {

            var mass = 0;
            var aa_obj = new AminoAcid();
            if(sequence) {
                for(var i = 0; i < sequence.length; i++) {
                    var aa = aa_obj.get(sequence.charAt(i));
                    mass += aa.avg;
                }
            }

            mass = _addTerminalModMass(mass, "n");
            mass = _addTerminalModMass(mass, "c");
            mass = _addResidueModMasses(mass, sequence.length, "n");
            // add N-terminal H
            mass = mass + Ion.MASS_H;
            // add C-terminal OH
            mass = mass + Ion.MASS_O + Ion.MASS_H;

            return mass;
        }

    this.getPotentialLosses = function _getPotentialLosses(sion)
    {
        var term = sion.term;
        var fragmentIndex = sion.fragmentIndex;
        if(term == 'n')
        {
            return nterm_totalLossOptions[fragmentIndex - 1];
        }
        if(term == 'c')
        {
            return cterm_totalLossOptions[sequence.length - fragmentIndex];
        }
    }

    this.customPotentialLosses = potentialLosses_custom;
    this.lorikeetPotentialLosses = potentialLosses_lorikeet;
    this.getLossForLabel = function _getLossForLabel(lossLabel)
    {
        var loss = potentialLosses_custom[lossLabel];
        if(loss) return loss;
        loss = potentialLosses_lorikeet[lossLabel];
        if(loss) return loss;
    }
    this.recalculateLossOptions = function _recalculate(selectedLossOptions, maxLossCount)
    {
        _calculateTotalLossOptions(selectedLossOptions, maxLossCount);
    }

    //----------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------

    function _initStaticMods()
    {
        if(staticModifications) {
            for(var i = 0; i < staticModifications.length; i += 1) {
                var mod = staticModifications[i];
                staticMods[mod.aa.code] = mod;
            }

            for(var i = 0; i < sequence.length; i += 1) {
                var mod = staticMods[sequence.charAt(i)];
                if(mod && mod.losses) {
                    for(var j = 0; j < mod.losses.length; j++) {
                        _addCustomLoss(i, mod.losses[j]);
                    }
                }
            }
        }
    }

    function _initVarMods()
    {
        if(varModifications) {
            for(var i = 0; i < varModifications.length; i += 1) {
                var mod = varModifications[i];
                if(mod) {
                    varMods[mod.position] = mod;

                    if(mod.losses) {
                        for(var j = 0; j < mod.losses.length; j++) {
                            var modIndex_0based = mod.position - 1; // mod.position is a 1-based index
                            _addCustomLoss(modIndex_0based, mod.losses[j]);
                        }
                    }
                }
            }
        }
    }

    function _addCustomLoss(index, loss)
    {
        // Example: {avgLossMass: "97.995", monoLossMass: "97.977", formula: "H3PO4"}
        if(loss && loss.avgLossMass && loss.avgLossMass > 0.0 && loss.monoLossMass && loss.monoLossMass > 0.0)
        {
            var neutralLoss = new NeutralLoss(loss.avgLossMass, loss.monoLossMass, loss.formula, loss.label);
            potentialLossesAtIndex[index].push(neutralLoss);
            potentialLosses_custom[neutralLoss.label()] = neutralLoss;
        }
    }

    function _initNeutralLosses()
    {
        for(var i = 0; i < sequence.length; i++) {
            potentialLossesAtIndex[i] = []; // potential neutral loss possibilities at an index in the peptide sequence
        }

        var ammoniaLoss = NeutralLoss.AmmoniaLoss();
        var waterLoss = NeutralLoss.WaterLoss();
        var phosphoLoss = NeutralLoss.PhosphoLoss();
        potentialLosses_lorikeet[ammoniaLoss.label()] = ammoniaLoss;
        potentialLosses_lorikeet[waterLoss.label()] = waterLoss;
        potentialLosses_lorikeet[phosphoLoss.label()] = phosphoLoss;

        for(var i = 0; i < sequence.length; i += 1)
        {
            var aa = sequence.charAt(i);
            if(aa == 'K' || aa == 'R' || aa == 'Q' || aa == 'N')
            {
                potentialLossesAtIndex[i].push(ammoniaLoss);
            }
            if (aa == 'S' || aa == 'T' || aa == 'E' || aa == 'D')
            {
                potentialLossesAtIndex[i].push(waterLoss);
            }
            if (aa == 'S' || aa == 'T' || aa == 'Y')
            {
                potentialLossesAtIndex[i].push(phosphoLoss);
            }
        }
    }

    function _calculateTotalLossOptions(selectedLosses, maxLossCount)
    {
        if(!maxLossCount)
            maxLossCount = 1;


        var selectedLossLabels = {};
        if(!selectedLosses)
        {
            for(var lossLabel in potentialLosses_custom)
            {
                selectedLossLabels[lossLabel] = true;
            }
        }
        else
        {
            for(var i = 0; i < selectedLosses.length; i += 1)
            {
                if(selectedLosses[i])
                {
                    selectedLossLabels[selectedLosses[i].label()] = true;
                }
            }
        }
        _initNeutralLossArrays(maxLossCount);
        _calculateTotalLossOptionsForTerm('n', selectedLossLabels, maxLossCount);
        _calculateTotalLossOptionsForTerm('c', selectedLossLabels, maxLossCount);

        if(debug)
        {
            _printPotentialLosses();
            _printNeutralLossCombinations();
        }
    }

    function _initNeutralLossArrays(maxLossCount)
    {
        for(var i = 0; i < sequence.length; i += 1)
        {
            nterm_totalLossOptions[i] = [];
            cterm_totalLossOptions[i] = [];
            for(var j = 0; j <= maxLossCount; j += 1)
            {
                nterm_totalLossOptions[i][j] = [];
                cterm_totalLossOptions[i][j] = [];
            }
        }
    }

    function _calculateTotalLossOptionsForTerm(term, selectedLossLabels, maxLossCount)
    {
        if(term == 'n')
        {
            for(var i = 0; i < sequence.length; i += 1) {
                _calculate(i, 0, 1, selectedLossLabels, nterm_totalLossOptions, maxLossCount);
            }
        }
        if(term == 'c')
        {
            var s = sequence.length - 1;
            for(var i = s; i >= 0; i -= 1) {
                _calculate(i, s, -1, selectedLossLabels, cterm_totalLossOptions, maxLossCount);
            }
        }
    }

    // totalLossOptionsArray[i] = lossOptions at index i in sequence.
    // totalLossOptionsArray[i][j] = LossCombinationList with j losses at index i in the sequence.
    // Each LossCombinationList contains an array of possible LossCombinations of j losses;
    function _calculate(i, s, incr, selectedLossLabels, totalLossOptionsArray, maxLossCount) {

        var lossOptions = potentialLossesAtIndex[i];
        var validLossOptions = [];
        for (var l = 0; l < lossOptions.length; l += 1) {  // these are the potential losses at index i in the sequence.
            var loss = lossOptions[l];
            if(!(loss.label() in selectedLossLabels))
                continue;
            else
                validLossOptions.push(loss);
        }

        if(validLossOptions.length == 0)
        {
            if (i != s)
            {
                totalLossOptionsArray[i] = totalLossOptionsArray[i - incr];
            }
            else
            {
                for(var j = 1; j <= maxLossCount; j += 1)
                {
                    totalLossOptionsArray[i][j] = new LossCombinationList(j);
                }
            }
            return;
        }

        var loss_options_at_i = [];
        totalLossOptionsArray[i] = loss_options_at_i;
        if(i == s)
        {
            for(j = 1; j <= maxLossCount; j += 1)
            {
                loss_options_at_i[j] = new LossCombinationList(j);
            }
        }
        else
        {
            var loss_options_before_i = totalLossOptionsArray[i - incr];
            for(j = 1; j <= maxLossCount; j += 1)
            {
                loss_options_at_i[j] = LossCombinationList.copyLossCombinationList(loss_options_before_i[j]);
            }
        }

        for(var l = 0; l < validLossOptions.length; l += 1)
        {
            var loss = validLossOptions[l];

            // iterate in reverse order so that we don't add the loss
            // in each iteration.
            for(var j = maxLossCount; j >= 1; j -= 1)
            {
                if(j == 1)
                {
                    var loss_combi_list_at_j = loss_options_at_i[j];
                    var loss_combi = new LossCombination();
                    loss_combi.addLoss(loss);
                    loss_combi_list_at_j.addLossCombination(loss_combi);
                }
                else
                {
                    var loss_combi_list_at_jminus1 = loss_options_at_i[j-1];
                    var loss_combi_list_at_j = LossCombinationList.copyLossCombinationList(loss_combi_list_at_jminus1);
                    loss_options_at_i[j] = loss_combi_list_at_j;
                    for(var k = 0; k < loss_combi_list_at_j.lossCombinationCount(); k += 1)
                    {
                        loss_combi_list_at_j.getLossCombination(k).addLoss(loss);
                    }
                }
            }
        }
    }

    function _getMassesWithLoss(lossOptions, mass) {
        var massesWithLoss = [];
        var j = 0;
        for (var i = 0; i < lossOptions.length; i += 1) {
            var loss = lossOptions[i];
            if (loss > 0.0) {
                var massWithLoss = mass - loss;
                massesWithLoss[j] = {mass:massWithLoss, loss:loss};
                j++;
            }
        }
        return massesWithLoss;
    }

    function _addResidueModMasses(seqMass, index, term) {

        var mass = seqMass;
        var slice = new _Slice(index, term);
        for( var i = slice.from; i < slice.to; i += 1) {
            // add any static modifications
            var mod = staticMods[sequence.charAt(i)];
            if(mod) {
                mass += mod.modMass;
            }
            // add any variable modifications
            mod = varMods[i+1]; // varMods index in the sequence is 1-based
            if(mod) {
                mass += mod.modMass;
            }
        }

        return mass;
    }

    function _getSeqMass(index, term, massType) {

        var mass = 0;
        var aa_obj = new AminoAcid();
        if(sequence) {
            var slice = new _Slice(index, term);
            for( var i = slice.from; i < slice.to; i += 1) {
                var aa = aa_obj.get(sequence[i]);
                mass += aa[massType];
            }
        }

        mass = _addTerminalModMass(mass, term);
        mass = _addResidueModMasses(mass, index, term);
        return mass;
    }

    function _addTerminalModMass(seqMass, term) {

        var mass = seqMass;
        // add any terminal modifications
        if(term == "n" && ntermMod)
            mass += ntermMod;
        if(term == "c" && ctermMod)
            mass += ctermMod;

        return mass;
    }

    function _Slice(index, term) {
        if(term == "n") {
            this.from = 0;
            this.to = index;
        }
        if(term == "c") {
            this.from = index;
            this.to = sequence.length;
        }
    }

    function _printNeutralLossCombinations() {

        for(var i = 1; i < sequence.length; i += 1)
        {
            var subseq = sequence.substring(0, i);
            var lossOpts = nterm_totalLossOptions[i];
            _log(subseq + " -- ");
            for(var j = 1; j < lossOpts.length; j += 1)
            {
                _log(" -- "+j);
                var lossesAt_j = lossOpts[j];
                var count = lossesAt_j.lossCombinationCount();
                for (var k = 0; k < count; k += 1)
                {
                    var lossCombo = lossesAt_j.getLossCombination(k);
                    _log("---- " + lossCombo.getLabels());
                }
            }
        }

        for(var i = sequence.length - 1; i >= 0; i -= 1)
        {
            var subseq = sequence.substring(i, sequence.length);
            var lossOpts = cterm_totalLossOptions[i];
            _log(subseq + " -- ");
            for(var j = 1; j < lossOpts.length; j += 1)
            {
                _log(" -- "+j);
                var lossesAt_j = lossOpts[j];
                var count = lossesAt_j.lossCombinationCount();
                for (var k = 0; k < count; k += 1)
                {
                    var lossCombo = lossesAt_j.getLossCombination(k);
                    _log("---- " + lossCombo.getLabels());
                }
            }
        }
    }

    function _printPotentialLosses() {
        for (var i = 0; i < potentialLossesAtIndex.length; i += 1) {
            var losses = potentialLossesAtIndex[i];
            if (losses.length > 0)
                _log("Potential Losses at " + i + " " + sequence.charAt(i));
            for (var j = 0; j < losses.length; j++) {
                var loss = losses[j];
                _log(loss.monoLossMass + ", " + loss.avgLossMass + ", " + loss.formula + ", " + loss.label());
            }
        }
    }

    function _log(message)
    {
        if(debug)
            console.log(message);
    }
}


//-----------------------------------------------------------------------------
// Modification
//-----------------------------------------------------------------------------
function Modification(aminoAcid, mass, losses) {
	this.aa = aminoAcid;
	this.modMass = mass;
    this.losses = losses;
}

function VariableModification(pos, mass, aminoAcid, losses) {
	this.position = parseInt(pos);
	this.aa = aminoAcid;
	this.modMass = mass;
    this.losses = losses;
}
//-----------------------------------------------------------------------------
// Neutral loss
//-----------------------------------------------------------------------------
function NeutralLoss(lossMassMono, lossMassAvg, formula, label)
{
    this.monoLossMass = lossMassMono;
    if(!this.monoLossMass)
        this.monoLossMass = 0.0;
    this.avgLossMass = lossMassAvg;
    if(!this.avgLossMass)
        this.avgLossMass = 0.0;
    this.formula = formula;
    this.userLabel = label;
    var _htmlLabel;

    this.longLabel = function _getLongLabel()
    {
        if(formula)
        {
            return formula + " ("+this.label()+")";
        }
        else
        {
            return this.label;
        }
    }

    this.label = function _getLabel()
    {
        if(this.userLabel) return this.userLabel;
        return "-"+Math.round(this.monoLossMass);
    }

    this.htmlLabel = function _getHtmlLabel()
    {
        if(_htmlLabel)
            return _htmlLabel;

        _htmlLabel = "";
        // html += H<sub>2</sub>O (<span style="font-weight: bold;">o</span>)
        if(formula)
        {
            for (var i = 0; i < formula.length; i++)
            {
                var charAt = formula.charAt(i);
                if(!isNaN(charAt))
                {
                    _htmlLabel += '<sub>'+charAt+'</sub>';
                }
                else
                {
                    _htmlLabel += charAt;
                }
            }
            _htmlLabel += " (";
        }
        _htmlLabel += '<span style="font-weight: bold;">'+this.label()+'</span>';
        if(formula)
        {
            _htmlLabel += ")";
        }

        return _htmlLabel;
    }
}

NeutralLoss.AmmoniaLoss = function _getAmmoniaLoss()
{
    return new NeutralLoss(Ion.AmmoniaLossMass_mono, Ion.AmmoniaLossMass_avg, "NH3", "*");
}
NeutralLoss.WaterLoss = function _getWaterLoss()
{
    return new NeutralLoss(Ion.WaterLossMass_mono, Ion.WaterLossMass_avg, "H2O", "o");
}
NeutralLoss.PhosphoLoss = function _getPhosphoLoss()
{
    return new NeutralLoss(Ion.PhosphoLossMass_mono, Ion.PhosphoLossMass_avg, "H3PO4", "p");
}

function LossCombination()
{
    this.losses = [];

    this.addLoss = function _addLoss(loss)
    {
        this.losses.push(loss);
    }

    this.getLabels = function _getLabels()
    {
        var label = "";
        for(var l = 0; l < this.losses.length; l += 1)
        {
            label += ", "+this.losses[l].label();
        }
        return label;
    }

    this.getLabel = function _getLabel()
    {
        if(this.losses.length == 1)
        {
            return " " + this.losses[0].label();
        }
        else
        {
            var lossMass = 0;
            for(var i = 0; i < this.losses.length; i += 1)
            {
                lossMass += this.losses[i].monoLossMass;
            }
            return " -" + Math.round(lossMass);
        }
    }

    this.getTotalLossMass = function _getTotalLossMass(massType)
    {
        var totalLoss = 0;
        for(var l = 0; l < this.losses.length; l += 1)
        {
            if(massType == 'mono')
            {
                totalLoss += this.losses[l].monoLossMass;
            }
            else if(massType == 'avg')
            {
                totalLoss += this.losses[l].avgLossMass;
            }
        }
        return totalLoss;
    }
}

LossCombination.copyLossCombination = function _copyLossCombination(originalLossCombination)
{
    var newLossCombination = new LossCombination();
    var losses = originalLossCombination.losses;
    for(var l = 0; l < losses.length; l += 1)
    {
        newLossCombination.addLoss(losses[l]);
    }
    return newLossCombination;
}

function LossCombinationList(numLosses)
{
    this.numLosses = numLosses;
    this.lossCombinations = [];

    this.lossCombinationCount = function _getLossCombinationCount()
    {
        return this.lossCombinations.length;
    }

    this.addLossCombination = function _addLossCombination(lossCombination)
    {
        var found = false;

        for(var l = 0; l < this.lossCombinationCount(); l += 1)
        {
            if(this.lossCombinations[l].getTotalLossMass('mono') == lossCombination.getTotalLossMass('mono'))
            {
                found = true;
            }
        }
        if(!found)
        {
            this.lossCombinations.push(lossCombination);
        }
    }

    this.getLossCombination = function _getLossCombinationAtIndex(index)
    {
        if(index >= 0 && index < this.lossCombinationCount())
        {
            return this.lossCombinations[index];
        }
    }
}
LossCombinationList.copyLossCombinationList = function _copyLossCombinationList(originalList)
{
    var newLossCombinationList = new LossCombinationList(originalList.numLosses);
    var count = originalList.lossCombinationCount();
    for(var l = 0; l < count; l += 1)
    {
        var originalLossCombination = originalList.getLossCombination(l);
        var lossCombination = LossCombination.copyLossCombination(originalLossCombination);
        newLossCombinationList.addLossCombination(lossCombination);
    }
    return newLossCombinationList;
}



