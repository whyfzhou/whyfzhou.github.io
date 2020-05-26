window.onload = init;
var selectedRows = [];

function init() {
  // bind events
  document.querySelector("#action-findss").onclick = findSteadyState;
  document.querySelector("#vload").onchange = vloadChanged;
  document.querySelector("#iload").onchange = iloadChanged;
  document.querySelector("#pload").onchange = ploadChanged;
  document.querySelector("#result-inverse").onclick = resultInverseClicked;
  document.querySelector("#result-remove").onclick = resultRemoveClicked;
  for (let i = 0; i < document.querySelector("#header").cells.length; ++i) {
    document.querySelector("#header").cells[i].onclick = showExplanation;
  }
  // document.querySelector("#measurementTable").onclick = tableClicked;
}

function showExplanation() {
  let explanation = this.innerHTML + ": " + this.getAttribute("data-tooltip");
  document.querySelector("#result-explanation").innerHTML = explanation;
  MathJax.Hub.Typeset(document.querySelector("#result-explanation"));
  // output = document.querySelector("#result-explanation");
  // MathJax.texReset();
  // let options = MathJax.getMetricsFor(output);
  // options.display = display.checked;
  // MathJax.tex2chtmlPromise(input, options)
  //   .then((node) => {
  //     output.appendChild(node);
  //     MathJax.startup.document.clear();
  //     MathJax.startup.document.updateDocument();
  //   })
  //   .catch((err) => {
  //     output.appendChild(document.createElement("pre")).appendChild(document.createTextNode(err.message));
  //   });
}

function refreshTable() {
  let table = document.querySelector("#rows");
  // let header = document.querySelector("#header");
  for (let row = 0; row < table.rows.length; ++row) {
    if (selectedRows.indexOf(row) === -1) {
      for (let col = 0; col < table.rows[row].cells.length; ++col) {
        table.rows[row].cells[col].classList.remove("table-selected");
      }
    } else {
      for (let col = 0; col < table.rows[row].cells.length; ++col) {
        table.rows[row].cells[col].classList.add("table-selected");
      }
    }
  }
}

function tableClicked() {
  let selected = this.rowIndex - 1;
  let j = selectedRows.indexOf(selected);
  if (j === -1) {
    selectedRows[selectedRows.length] = selected;
  } else {
    selectedRows.splice(j, 1);
  }
  // console.log(selectedRows);
  refreshTable();
}

function resultInverseClicked() {
  let table = document.querySelector("#rows");
  for (let i = 0; i < table.rows.length; ++i) {
    let j = selectedRows.indexOf(i);
    if (j === -1) {
      selectedRows[selectedRows.length] = i;
    } else {
      selectedRows.splice(j, 1);
    }
  }
  // console.log(selectedRows);
  refreshTable();
}

function resultRemoveClicked() {
  selectedRows.sort();
  let table = document.querySelector("#rows");
  for (let i = selectedRows.length - 1; i >= 0; --i) {
    table.deleteRow(selectedRows[i]);
  }
  if (table.rows.length === 0) {
    document.querySelector("#header").setAttribute("style", "visibility: hidden;");
    for (let el of document.querySelectorAll(".result")) {
      el.setAttribute("style", "visibility: hidden;");
    }
  }
  selectedRows = [];
}

function vloadChanged() {
  let vload = document.querySelector("#vload").valueAsNumber;
  let iload = document.querySelector("#iload").valueAsNumber;
  let pload = vload * iload;
  document.querySelector("#pload").value = pload;
}

function iloadChanged() {
  let vload = document.querySelector("#vload").valueAsNumber;
  let iload = document.querySelector("#iload").valueAsNumber;
  let pload = vload * iload;
  document.querySelector("#pload").value = pload;
}

function ploadChanged() {
  let vload = document.querySelector("#vload").valueAsNumber;
  let pload = document.querySelector("#pload").valueAsNumber;
  let iload = pload / vload;
  document.querySelector("#iload").value = iload;
}

function findSteadyState() {
  let lr = document.querySelector("#lr").valueAsNumber * 1e-6;
  let lm = document.querySelector("#lm").valueAsNumber * 1e-6;
  let cr = document.querySelector("#cr").valueAsNumber * 1e-9;
  let nps = document.querySelector("#nps").valueAsNumber;
  let t12 = document.querySelector("#t12").valueAsNumber * 1e-9;
  let vbus = document.querySelector("#vbus").valueAsNumber;
  let vload = document.querySelector("#vload").valueAsNumber;
  let pload = document.querySelector("#pload").valueAsNumber;
  let vout = vload * nps;

  let result = steady_state_pout(pload, lr, lm, cr, vbus, vout, t12);

  let data = sample(result, 2000);
  let measurement = evaluateSampled(result, data);
  fillTable(measurement);

  drawSingleOP(data);
}

function fillTable(measurement) {
  for (let resultAction of document.querySelectorAll(".result-action")) {
    resultAction.setAttribute("style", "visibility: visable");
  }
  for (let resultAction of document.querySelectorAll(".result")) {
    resultAction.setAttribute("style", "visibility: visable");
  }
  document.querySelector("#header").setAttribute("style", "visibility: visible;");
  let nps = document.querySelector("#nps").valueAsNumber;
  let row = [
    document.querySelector("#lr").valueAsNumber,
    document.querySelector("#lm").valueAsNumber,
    document.querySelector("#cr").valueAsNumber,
    document.querySelector("#nps").valueAsNumber,
    document.querySelector("#t12").valueAsNumber,
    document.querySelector("#vbus").valueAsNumber,
    document.querySelector("#vload").valueAsNumber,
    document.querySelector("#iload").valueAsNumber,
    document.querySelector("#pload").valueAsNumber,
    measurement.fsw / 1e3,
    measurement.dutyHighSide * 100,
    measurement.dutyDiode * 100,
    measurement.iRMSPri,
    measurement.vMeanRes,
    measurement.vACRMSRes,
    measurement.vRMSPri,
    measurement.iRMSOut * nps,
    measurement.iPPeakOut * nps,
    measurement.psiPri / 1e-3,
    measurement.qLowHigh / 1e-9,
  ];

  let newRow = document.createElement("tr");
  let header = document.querySelector("#header");
  let i = 0;
  for (let x of row) {
    let newCell = document.createElement("td");
    newCell.innerHTML = x.toLocaleString(undefined, {
      maximumFractionDigits: 2,
    });
    if (header.cells[i].className === "table-header-parameter") {
      newCell.classList.add("table-parameter");
    } else if (header.cells[i].className === "table-header-operatingpoint") {
      newCell.classList.add("table-operatingpoint");
    } else {
      newCell.classList.add("table-performance");
    }
    newRow.append(newCell);
    i++;
  }
  newRow.onclick = tableClicked;
  document.querySelector("#rows").appendChild(newRow);
}
