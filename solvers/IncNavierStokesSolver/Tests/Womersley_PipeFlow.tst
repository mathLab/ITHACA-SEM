<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Womersley B.C. Test</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Womersley_PipeFlow.xml</parameters>
    <files>
        <file description="Session File">Womersley_PipeFlow.xml</file>
        <file description="Session File">fourier_coef.txt</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">2.26003e-3</value>
            <value variable="v" tolerance="1e-8">1.99059e-3</value>
	        <value variable="w" tolerance="1e-8">3.60492e-2</value>
	        <value variable="p" tolerance="1e-8">9.74819e-2</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">1.48938e-2</value>
            <value variable="v" tolerance="1e-8">1.07501e-2</value>
            <value variable="w" tolerance="1e-8">1.33385e-1</value>
            <value variable="p" tolerance="1e-8">6.80866e-1</value>
        </metric>
    </metrics>
</test>

