const form = document.querySelector('#analysis-form');
const statusNode = document.querySelector('#status');
const results = document.querySelector('#results');
const cards = document.querySelector('#cards');

fetch('/api/v1/models')
  .then((response) => response.json())
  .then((models) => {
    const status = models.bundle.evaluation_status || 'evaluation status unavailable';
    document.querySelector('#model-status').textContent =
      `${models.bundle.split_genomes || 'Unknown'} grouped genomes · ${status}`;
  })
  .catch(() => {
    document.querySelector('#model-status').textContent = 'Model provenance unavailable';
  });

form.addEventListener('submit', async (event) => {
  event.preventDefault();
  results.hidden = true;
  statusNode.textContent = 'Running assembly QC, M1, M2, and model inference…';
  const data = new FormData();
  data.append('fasta', document.querySelector('#fasta').files[0]);
  try {
    const response = await fetch('/api/v1/analyses', { method: 'POST', body: data });
    const report = await response.json();
    if (!response.ok) throw new Error(report.detail || 'Analysis failed');
    cards.replaceChildren(...report.decisions.map((decision) => {
      const node = document.createElement('article');
      node.className = 'card';
      node.innerHTML = `<h3>${decision.antibiotic}</h3>
        <p class="call">${decision.call.replaceAll('_', ' ')}</p>
        <p>${(decision.confidence * 100).toFixed(1)}% confidence</p>
        <p>${decision.evidence_category.replaceAll('_', ' ')}</p>
        <p>Target: ${decision.target_status}</p>`;
      return node;
    }));
    document.querySelector('#evidence').textContent = JSON.stringify(report.workflows, null, 2);
    results.hidden = false;
    statusNode.textContent = `Analysis ${report.analysis_id} complete.`;
  } catch (error) {
    statusNode.textContent = error.message;
  }
});
