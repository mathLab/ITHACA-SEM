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
            <value variable="u" tolerance="1e-12">5.47382e-08</value>
            <value variable="v" tolerance="1e-12">5.02179e-08</value>
            <value variable="p" tolerance="1e-12">2.33967e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">4.01658e-08</value>
            <value variable="v" tolerance="1e-12">3.98244e-08</value>
            <value variable="p" tolerance="1e-12">1.60286e-06</value>
        </metric>
    </metrics>
</test>