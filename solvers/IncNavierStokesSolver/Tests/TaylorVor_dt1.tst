<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Convergence: Taylor Vortex IMEXOrder2 dt=0.01</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>TaylorVor_dt1.xml</parameters>
    <files>
        <file description="Session File">TaylorVor_dt1.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">5.95129e-06</value>
            <value variable="v" tolerance="1e-12">4.9951e-06</value>
            <value variable="p" tolerance="1e-12">0.000232625</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">4.1544e-06</value>
            <value variable="v" tolerance="1e-12">4.51559e-06</value>
            <value variable="p" tolerance="1e-12">0.00016244</value>
        </metric>
    </metrics>
</test>