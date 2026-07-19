from streamlit.testing.v1 import AppTest


def test_antibiotic_selector_survives_reruns() -> None:
    app = AppTest.from_file("sandbox/app.py", default_timeout=30).run()
    assert not app.exception

    observation = next(widget for widget in app.selectbox if widget.label == "Antibiotic")
    diagnostics = next(
        widget for widget in app.selectbox if widget.label == "Model diagnostics antibiotic"
    )
    for value in observation.options:
        observation.select(value)
        diagnostics.select(value)
        app.run()
        assert not app.exception
