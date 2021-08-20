var cols = ['CID', 'SMILES', 'AllAssays_Active', 'AllAssays_Total',
 'LuciferaseAssays_Active', 'LuciferaseAssays_Total', 'BetaLactamaseAssays_Active',
 'BetaLactamaseAssays_Total', 'FluorescenceAssays_Active', 'FluorescenceAssays_Total', 
 'PAINS Alert', 'Tc'];
 
function createHighlightImg(highlights, smarts_or_smile, smile){
	var highlightSrc = "/highlighted-image/" + smarts_or_smile + "/";
	for(var i = 0; i < highlights.length-1; i++){
			highlightSrc += highlights[i] + "-";
	}
	var re = new RegExp('/', 'g');
	smile = smile.replace(re, "FOWARDSLASH");
	highlightSrc += highlights[highlights.length-1] + "_" + smile + ".png";
	return escape(highlightSrc);
}
function clearDisplay(){
	$("#noPAINS").text("");
	$("#similarCSVLink").text("");
	$("#painsCSVLink").text("");
	$("#noPAINS").text("");
	$("#pains_table").find("tr:gt(0)").remove();
	$("#img_table").find("tr:gt(1)").remove();
//	$("#same_pains_A_table").find("tr:gt(0)").remove();
	$("#other_pains_table").find("tr:gt(0)").remove();
	$("#non_pains_table").find("tr:gt(0)").remove();
	$("#pains_container").empty();
	$("#painsCSVLink").removeClass();
//	$("#same_pains_A_csv").removeClass();
	$("#other_pains_csv").removeClass();
	$("#non_pains_csv").removeClass();
	$("#all_csv").removeClass();
	
	$("#painsCSVLink").text("");
//	$("#same_pains_A_csv").text("");
	$("#other_pains_csv").text("");
	$("#non_pains_csv").text("");
	$("#all_csv").text("");
}

const col = ['Overall', 'NN', 'Tc', 'Name', 'Exp_BioAct'];


function makeRow(dict){
	var row = $("<tr/>");
	row.prop("colspan","4");
	//row.append('<td> <a href="https://pubchem.ncbi.nlm.nih.gov/compound/'+dict['cid']+'">'+dict['cid']+'</a></td>');
	//var img = "/image/smile/" + dict['smile'] + '.png';
	//if(pains){
	//	img = "/highlighted-image/smile/" + dict['pains_highlights'] + '_' + dict['smile'] + '.png';
	//}
	//row.append('<td> <img src="'+escape(img)+'"></td>');
	
	row.append('<td>'+dict[col[0]]+'</td>');
	row.append('<td>'+dict[col[1]]+'</td>');
	var roundedTc = Math.round(dict["tc"] * 100) / 100;
	row.append('<td>'+roundedTc +'</td>');
	row.append('<td>'+dict[col[3]]+'</td>');
	row.append('<td>'+dict[col[4]]+'</td>');

	return row;
}

function makeCSVRow(dict, pains){
	csvRow = dict['cid'] + ',' + dict['smile'] + ',' +  (Math.round(dict["tc"] * 100) / 100);
	for(var i = 0; i < col.length; i++){
		 csvRow+= ',' + dict[col[i]];
	}
	if(pains){
		csvRow+= ',' + dict['pains'];
	}
	csvRow += "\r\n";
	return csvRow;
 
}

function buildPainsTable(id){
	var t = '<table align="center" id="'+ id +'_table" class="table table-bordered">';
	t += '<tr colspan="4">';
	t += '<th class="col-md-1">CID</th>';
	t += '<th class="col-md-1">Substructure</th>';
	t += '<th class="col-md-1">Tc</th>';
	t += '<th class="col-md-3">Analysis</th>';
	t += '<th class="col-md-1">All Assays Active</th>';
	t += '<th class="col-md-1">All Assays Total</th>';
	t += '<th class="col-md-1">Luciferase Assays Active</th>';
	t += '<th class="col-md-1">Luciferase Assays Total</th>';
	t += '<th class="col-md-1">BetaLactamase Assays Active</th>';
	t += '<th class="col-md-1">BetaLactamase Assays Total</th>';
	t += '<th  class="col-md-.5">Fluorescence Assays Active</th>';
	t += '<th  class="col-md-.5">Fluorescence Assays Total</th>';
	t += '<th  class="col-md-1">PAINS Alert</th>';
	t += '</tr>';
	t += '</table>';
    t += '<a id="'+id+'_csv"></a>';
	$("#pains_container").append(t);
}
var storedJson = null;

function martinsOutputToTable(smile_json){
	clearDisplay();
	if(smile_json=="Error, poorly formed smile string."){
		$("#noPAINS").text("Error, poorly formed smile string.");
		return;
	}
	storedJson = smile_json;
	//var most = smile_json["pains_most"];
	var maccs_nearest = smile_json["maccs"];
	var isida_nearest = smile_json["isida"];
	var dragon_nearest = smile_json["dragon"];

	//var pains = Object.keys(highlightsForQuery);
	
	for(var i = 0; i < maccs_nearest.length; i++){
		var dict = maccs_nearest[i];	
		if(typeof dict == 'object'){
			var row = makeRow(dict);
			$("#maccs_table").append(row);
		}
	}


	//var row = makeRow(smile_json["non_pains_most"],false);
	//$("#non_pains_table").append(row);
	
	var other = smile_json["non_pains_nearest"];
	for(var i = 0; i < other.length; i++){
		var dict = other[i];	
		if(typeof dict == 'object'){
			var row = makeRow(dict,false);
			$("#non_pains_table").append(row);
		}
	}


	for(var i = 0; i < painsAlerts.length; i++){ //appends painsAlerts to pains table
		var painsArray = painsAlerts[i];	
		if(typeof painsArray == 'object'){
			var row = $("<tr/>");
			row.append('<td>'+painsArray[0]+'</td>');
			var encoded = encodeURIComponent(painsArray[0])
			row.append('<td> <img src="/image/smarts/'+encoded+'.png"></td>')
			for( var j = 1; j < painsArray.length; j++){
				row.append('<td>'+painsArray[j]+'</td>');	
			}
			$("#pains_table").append(row);
		}
	}

	similarCompoundsCSV();
	painsCSV();
}


function checkPAINS() {
  smile = $('#jme_output').val()
  req = JSON.stringify({"smiles": smile});
  $.ajax({
    url:"_PAINS",
    type:"POST",
    data:req,
    contentType:"application/json; charset=utf-8",
    success: function(response, textStatus, jqXHR) {
		martinsOutputToTable(response);
    },
    error: function(jqXHR, textStatus, errorThrown){
      document.getElementById("script_output").innerHTML = 'Poorly formatted Smile String';
   }

  });

}

function similarCompoundsCSV() {
	var row_data = ['CID', 'SMILES', 'Tc', 'AllAssays_Active', 'AllAssays_Total',
 	'LuciferaseAssays_Active', 'LuciferaseAssays_Total', 'BetaLactamaseAssays_Active',
 	'BetaLactamaseAssays_Total', 'FluorescenceAssays_Active', 'FluorescenceAssays_Total', 
 	'PAINS Alert'];

	
	var allCsvContent = "data:text/csv;charset=utf-8,";
	allCsvContent += row_data.join(",");
	allCsvContent += "\r\n";//add headers

	var sameCsvContent = "data:text/csv;charset=utf-8,";
	sameCsvContent += row_data.join(",");
	sameCsvContent += "\r\n";//add headers

	
	var otherCsvContent = "data:text/csv;charset=utf-8,";
	otherCsvContent += row_data.join(",");
	otherCsvContent += "\r\n";//add headers
	
	var nonCsvContent = "data:text/csv;charset=utf-8,";
	row_data.pop();//get rid of PAINS alert
	nonCsvContent += row_data.join(",");
	nonCsvContent += "\r\n";//add headers

	
	//we need to generate the CSV links for the pains matches dynamically.
	//This is done when we loop over all the different PAINS matches to make images.
	for(var i = 0; i < storedJson["pains_nearest"].length; i++){
		var sameCsvContent = "data:text/csv;charset=utf-8,";
		sameCsvContent += row_data.join(",");
		sameCsvContent += "\r\n";//add headers
		for(var j = 0; j < storedJson["pains_nearest"][i].length; j++){
			allCsvContent += makeCSVRow(storedJson["pains_nearest"][i][j], true);
			sameCsvContent += makeCSVRow(storedJson["pains_nearest"][i][j], true);
		}
		var sameEncodedUri = encodeURI(sameCsvContent);
		var csv_name = "pains_" + i +"_csv";
		$(document.getElementById(csv_name)).prop("href",sameEncodedUri);
		$(document.getElementById(csv_name)).prop("download","similarCompoundWithPAINS.csv");
		$(document.getElementById(csv_name)).text("Download Similar Compounds with matching PAINS Alerts as CSV");
		$(document.getElementById(csv_name)).addClass("btn btn-primary btn-lg");
	}
	
	
	for(var i = 0; i < storedJson["other_pains_nearest"].length; i++){
		allCsvContent += makeCSVRow(storedJson["other_pains_nearest"][i], true);
		otherCsvContent += makeCSVRow(storedJson["other_pains_nearest"][i], true);
	}
	
	for(var i = 0; i < storedJson["non_pains_nearest"].length; i++){
		allCsvContent += makeCSVRow(storedJson["non_pains_nearest"][i], false); 
		nonCsvContent += makeCSVRow(storedJson["non_pains_nearest"][i], false);
	}
	allCsvContent += '\r\n';
	allCsvContent += painsCSVText();

	var allEncodedUri = encodeURI(allCsvContent);
	
	var otherEncodedUri = encodeURI(otherCsvContent);
	var nonEncodedUri = encodeURI(nonCsvContent);
	



	$("#other_pains_csv").prop("href",otherEncodedUri);
	$("#other_pains_csv").prop("download","similarCompoundWithDifferentPAINS.csv");
	$("#other_pains_csv").text("Download Similar Compounds without matching Pains as CSV");
	$("#other_pains_csv").addClass("btn btn-primary btn-lg");
	

	$("#non_pains_csv").prop("href",nonEncodedUri);
	$("#non_pains_csv").prop("download","similarCompoundWithoutPAINS.csv");
	$("#non_pains_csv").text("Download Similar Compounds with No Pains as CSV");
	$("#non_pains_csv").addClass("btn btn-primary btn-lg");

	$("#all_csv").prop("href",allEncodedUri);
	$("#all_csv").prop("download","similarCompoundsAndPainsAlerts.csv");
	$("#all_csv").text("Download Similar Compounds and PAINS Alerts as CSV");
	$("#all_csv").addClass("btn btn-primary btn-lg");

// Citation:
// https://stackoverflow.com/questions/14964035/how-to-export-javascript-array-info-to-csv-on-client-side}
}


function painsCSVText() {
	var row_data = ['Smarts','PAINS_Alert','N_PubChem','N_DCM', 'Luciferase', 'BetaLactamase','Fluorescence','All_Assays'];

	csvContent = row_data.join(',');
	csvContent += '\r\n';

	pains = storedJson["flagged_alerts"];
	for(var i = 0; i < pains.length; i++){
		csvContent += '"' + pains[i][0] + '",';
		csvContent += pains[i].slice(1).join(',');
		csvContent += '\r\n';
	}
	return csvContent; 
}

function painsCSV(){
	var csvContent = "data:text/csv;charset=utf-8,";
	csvContent += painsCSVText();
	
	var encodedUri = encodeURI(csvContent);
	$("#painsCSVLink").prop("href",encodedUri);
	$("#painsCSVLink").prop("download","PAINSAlerts.csv");
	$("#painsCSVLink").text("Download PAINS Alerts as CSV");
	$("#painsCSVLink").addClass("btn btn-primary btn-lg");
	
	// Citation:
	// https://stackoverflow.com/questions/14964035/how-to-export-javascript-array-info-to-csv-on-client-side}
}
