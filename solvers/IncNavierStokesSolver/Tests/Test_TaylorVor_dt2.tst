<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Convergence: Taylor Vortex IMEXOrder2 dt=0.001</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Test_TaylorVor_dt2.xml</parameters>
    <files>
        <file description="Session File">Test_TaylorVor_dt2.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">4.3947e-08</value>
            <value variable="v" tolerance="1e-12">1.92721e-07</value>
            <value variable="p" tolerance="1e-12">2.54032e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">3.3274e-08</value>
            <value variable="v" tolerance="1e-12">1.04068e-07</value>
            <value variable="p" tolerance="1e-12">5.96329e-06</value>
        </metric>
    </metrics>
</test>