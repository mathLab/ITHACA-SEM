<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D equilateral triangle flow, Tetrahedral elements, P=5</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_equitri.xml</parameters>
    <files>
        <file description="Session File">Tet_equitri.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-06">3.7863e-11</value>
            <value variable="v" tolerance="1e-06">1.26062e-10</value>
            <value variable="w" tolerance="1e-06">7.94176e-09</value>
	    <value variable="p" tolerance="1e-06">8.54695e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-06">7.13781e-10</value>
            <value variable="v" tolerance="1e-06">1.81283e-09</value>
            <value variable="w" tolerance="1e-06">2.59194e-08</value>
	    <value variable="p" tolerance="1e-06">3.84038e-07</value>
        </metric>
    </metrics>
</test>
