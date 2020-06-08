window.onload = init;
window.app = [];
window.app.selectedRows = [];
window.app.operatingPointOrder = [0, 1, 2, 3, 4];
window.app.operatingPointNumber = 1;
window.app.maximumOperatingPoint = 5;
window.app.resultIndices = [];

// ----------------------------------------------------------------------------

function init() {
  // bind events
  document.querySelector("#action-findss").onclick = findSteadyState0;
  document.querySelector("#action-findssa").onclick = findSteadyStateAll;
  document.querySelector("#vload").onchange = vloadChanged;
  document.querySelector("#iload").onchange = iloadChanged;
  document.querySelector("#pload").onchange = ploadChanged;
  document.querySelector("#result-inverse").onclick = resultInverseClicked;
  document.querySelector("#result-remove").onclick = resultRemoveClicked;
  document.querySelector("#action-addremoveop").onclick = addOperatingPoint;
  for (
    let i = 0;
    i < document.querySelector("#result-table-header").cells.length;
    ++i
  ) {
    document.querySelector("#result-table-header").cells[
      i
    ].onclick = showExplanation;
  }
}

// ----------------------------------------------------------------------------

function addOperatingPoint() {
  if (window.app.operatingPointNumber >= window.app.maximumOperatingPoint) {
    return;
  }
  let newOpIndex =
    window.app.operatingPointOrder[window.app.operatingPointNumber];
  let opTemplate = document.querySelector("#operating-point-template");
  let newOp = opTemplate.content.cloneNode(true);
  for (let n of newOp.querySelectorAll("*")) {
    let nid = n.getAttribute("id");
    if (nid) {
      n.setAttribute("id", `${nid}${newOpIndex}`);
    }
  }
  for (let n of newOp.querySelectorAll("label.unit")) {
    let nfor = n.getAttribute("for");
    n.setAttribute("for", `${nfor}${newOpIndex}`);
  }
  newOp.querySelector("div.op").setAttribute("id", `op${newOpIndex}`);
  let opContainer = document.querySelector("#operating-points-container");
  opContainer.appendChild(newOp);
  window.app.operatingPointNumber++;
  document.querySelector(
    `#action-addremoveop${newOpIndex}`
  ).onclick = removeOperatingPoint;
  document.querySelector(`#vload${newOpIndex}`).onclick = vloadChanged;
  document.querySelector(`#iload${newOpIndex}`).onclick = iloadChanged;
  document.querySelector(`#pload${newOpIndex}`).onclick = ploadChanged;
  document.querySelector(`#vload${newOpIndex}`).value = document.querySelector(
    "#vload"
  ).valueAsNumber;
  document.querySelector(`#iload${newOpIndex}`).value = document.querySelector(
    "#iload"
  ).valueAsNumber;
  document.querySelector(`#pload${newOpIndex}`).value = document.querySelector(
    "#pload"
  ).valueAsNumber;
}

function removeOperatingPoint() {
  let opIndex = parseInt(this.id.slice(-1)[0]);
  let opToDelete = document.querySelector(`#op${opIndex}`);
  opToDelete.parentNode.removeChild(opToDelete);
  window.app.operatingPointOrder.splice(
    window.app.operatingPointOrder.indexOf(opIndex),
    1
  );
  window.app.operatingPointOrder.push(opIndex);
  window.app.operatingPointNumber--;
}

function vloadChanged() {
  let opIndex = parseInt(this.id.slice(-1)[0]);
  if (Number.isNaN(opIndex)) {
    opIndex = "";
  }
  let vload = document.querySelector(`#vload${opIndex}`).valueAsNumber;
  let iload = document.querySelector(`#iload${opIndex}`).valueAsNumber;
  let pload = vload * iload;
  document.querySelector(`#pload${opIndex}`).value = pload;
}

function iloadChanged() {
  let opIndex = parseInt(this.id.slice(-1)[0]);
  if (Number.isNaN(opIndex)) {
    opIndex = "";
  }
  let vload = document.querySelector(`#vload${opIndex}`).valueAsNumber;
  let iload = document.querySelector(`#iload${opIndex}`).valueAsNumber;
  let pload = vload * iload;
  document.querySelector(`#pload${opIndex}`).value = pload;
}

function ploadChanged() {
  let opIndex = parseInt(this.id.slice(-1)[0]);
  if (Number.isNaN(opIndex)) {
    opIndex = "";
  }
  let vload = document.querySelector(`#vload${opIndex}`).valueAsNumber;
  let pload = document.querySelector(`#pload${opIndex}`).valueAsNumber;
  let iload = pload / vload;
  document.querySelector(`#iload${opIndex}`).value = iload;
}

// ----------------------------------------------------------------------------

function findSteadyState0() {
  findSteadyState(0);
}

function findSteadyStateAll() {
  for (let i = 0; i < window.app.operatingPointNumber; ++i) {
    findSteadyState(window.app.operatingPointOrder[i]);
  }
}

function findSteadyState(n) {
  let opIndex = !n || n === 0 ? "" : n;
  let lr = document.querySelector("#lr").valueAsNumber * 1e-6;
  let lm = document.querySelector("#lm").valueAsNumber * 1e-6;
  let cr = document.querySelector("#cr").valueAsNumber * 1e-9;
  let nps = document.querySelector("#nps").valueAsNumber;
  let t12 = document.querySelector("#t12").valueAsNumber * 1e-9;
  let vbus = document.querySelector(`#vbus${opIndex}`).valueAsNumber;
  let vload = document.querySelector(`#vload${opIndex}`).valueAsNumber;
  let pload = document.querySelector(`#pload${opIndex}`).valueAsNumber;
  let vout = vload * nps;

  let currentResultId;
  if (window.app.resultIndices.length === 0) {
    currentResultId = 1;
  } else {
    currentResultId = window.app.resultIndices.slice(-1)[0] + 1;
  }
  window.app.resultIndices.push(currentResultId);

  let sampledData;
  let measurement;
  try {
    let steadyState = steady_state_pout(pload, lr, lm, cr, vbus, vout, t12);
    sampledData = sample(steadyState, 2000);
    measurement = evaluateSampled(steadyState, sampledData);
  } catch (e) {
    let pmax = e.value;
    sampledData = pmax;
    measurement = pmax;
  } finally {
    fillTable(measurement, opIndex);
    drawSingleOP(sampledData);
  }
}

// ----------------------------------------------------------------------------

function setExplanation(s) {
  document
    .querySelector("#result-explanation")
    .classList.remove("result-explanation-default");
  document
    .querySelector("#result-explanation")
    .classList.add("result-explanation-clicked");
  document.querySelector("#result-explanation").innerHTML = s;
}

function resetExplanation() {
  document
    .querySelector("#result-explanation")
    .classList.add("result-explanation-default");
  document
    .querySelector("#result-explanation")
    .classList.remove("result-explanation-clicked");
  document.querySelector("#result-explanation").innerHTML =
    "Click head titles for more information.";
}

function showExplanation() {
  let explanation = this.innerHTML + ": " + this.getAttribute("data-tooltip");
  setExplanation(explanation);
  MathJax.Hub.Typeset(document.querySelector("#result-explanation"));
}

// ----------------------------------------------------------------------------

function fillTable(measurement, n) {
  let opIndex = !n || n === 0 ? "" : n;
  for (let resultAction of document.querySelectorAll(".result-action")) {
    resultAction.setAttribute("style", "display: inline-block;");
  }
  for (let resultAction of document.querySelectorAll(".result")) {
    resultAction.setAttribute("style", "display: inline;");
  }
  document
    .querySelector("#result-table-header")
    .setAttribute("style", "display: table-row;");
  let nps = document.querySelector("#nps").valueAsNumber;
  let row;
  if (!(typeof measurement === "number")) {
    row = [
      window.app.resultIndices.slice(-1)[0],
      document.querySelector("#lr").valueAsNumber,
      document.querySelector("#lm").valueAsNumber,
      document.querySelector("#cr").valueAsNumber,
      document.querySelector("#nps").valueAsNumber,
      document.querySelector("#t12").valueAsNumber,
      document.querySelector(`#vbus${opIndex}`).valueAsNumber,
      document.querySelector(`#vload${opIndex}`).valueAsNumber,
      document.querySelector(`#iload${opIndex}`).valueAsNumber,
      document.querySelector(`#pload${opIndex}`).valueAsNumber,
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
  } else {
    row = [
      window.app.resultIndices.slice(-1)[0],
      document.querySelector("#lr").valueAsNumber,
      document.querySelector("#lm").valueAsNumber,
      document.querySelector("#cr").valueAsNumber,
      document.querySelector("#nps").valueAsNumber,
      document.querySelector("#t12").valueAsNumber,
      document.querySelector(`#vbus${opIndex}`).valueAsNumber,
      document.querySelector(`#vload${opIndex}`).valueAsNumber,
      document.querySelector(`#iload${opIndex}`).valueAsNumber,
      document.querySelector(`#pload${opIndex}`).valueAsNumber,
    ];
  }

  let newRow = document.createElement("tr");
  let header = document.querySelector("#result-table-header");
  let i = 0;
  for (let x of row) {
    let newCell = document.createElement("td");
    if (i === 0) {
      let a = document.createElement("a");
      a.href = `#fig${window.app.resultIndices.slice(-1)[0]}`;
      a.innerHTML = x.toLocaleString(undefined, {
        maximumFractionDigits: 2,
      });
      newCell.appendChild(a);
    } else {
      newCell.innerHTML = x.toLocaleString(undefined, {
        maximumFractionDigits: 2,
      });
    }
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

  if (typeof measurement === "number") {
    let errorValue = measurement;
    let newCell = document.createElement("td");
    newCell.classList.add("result-error-message");
    newCell.setAttribute("colspan", "11");
    newCell.innerText = `${errorValue}`;
    newRow.append(newCell);
  }

  newRow.onclick = tableClicked;
  document.querySelector("#result-table-rows").appendChild(newRow);
}

function refreshTable() {
  let table = document.querySelector("#result-table-rows");
  for (let row = 0; row < table.rows.length; ++row) {
    let resultId = parseInt(table.rows[row].cells[0].innerText);
    let caption = document.querySelector(`#fig${resultId} p`);
    if (window.app.selectedRows.indexOf(row) === -1) {
      for (let col = 0; col < table.rows[row].cells.length; ++col) {
        table.rows[row].cells[col].classList.remove("result-selected");
        caption.classList.remove("result-selected");
      }
    } else {
      for (let col = 0; col < table.rows[row].cells.length; ++col) {
        table.rows[row].cells[col].classList.add("result-selected");
        caption.classList.add("result-selected");
      }
    }
  }
}

function tableClicked(e) {
  if (e.target.nodeName !== "A") {
    let selected = this.rowIndex - 1;
    let j = window.app.selectedRows.indexOf(selected);
    if (j === -1) {
      window.app.selectedRows[window.app.selectedRows.length] = selected;
    } else {
      window.app.selectedRows.splice(j, 1);
    }
    refreshTable();
  }
}

function resultInverseClicked() {
  let table = document.querySelector("#result-table-rows");
  for (let i = 0; i < table.rows.length; ++i) {
    let j = window.app.selectedRows.indexOf(i);
    if (j === -1) {
      window.app.selectedRows.push(i);
    } else {
      window.app.selectedRows.splice(j, 1);
    }
  }
  refreshTable();
}

function resultRemoveClicked() {
  window.app.selectedRows.sort((a, b) => a - b);
  let table = document.querySelector("#result-table-rows");
  for (let i = window.app.selectedRows.length - 1; i >= 0; --i) {
    let idToRemove = parseInt(
      table.rows[window.app.selectedRows[i]].cells[0].textContent
    );
    let figToRemove = document.querySelector(`#fig${idToRemove}`);
    figToRemove.parentNode.removeChild(figToRemove);
    table.deleteRow(window.app.selectedRows[i]);
  }
  if (table.rows.length === 0) {
    document
      .querySelector("#result-table-header")
      .setAttribute("style", "display: none;");
    for (let el of document.querySelectorAll(".result")) {
      el.setAttribute("style", "display: none;");
    }
    for (let resultAction of document.querySelectorAll(".result-action")) {
      resultAction.setAttribute("style", "display: none;");
    }
    for (let resultAction of document.querySelectorAll(".result")) {
      resultAction.setAttribute("style", "display: none;");
    }
    resetExplanation();
  }
  window.app.selectedRows = [];
}
