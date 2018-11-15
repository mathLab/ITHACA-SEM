<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Convergence: Taylor Vortex IMEXOrder2 dt=0.001</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>TaylorVor_dt2.xml</parameters>
    <files>
        <file description="Session File">TaylorVor_dt2.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">5.48624e-08</value>
            <value variable="v" tolerance="1e-12">5.00282e-08</value>
            <value variable="p" tolerance="5e-11">2.33974e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">4.02333e-08</value>
            <value variable="v" tolerance="1e-12">3.97501e-08</value>
            <value variable="p" tolerance="2e-11">1.60287e-06</value>
        </metric>
    </metrics>
</test>
