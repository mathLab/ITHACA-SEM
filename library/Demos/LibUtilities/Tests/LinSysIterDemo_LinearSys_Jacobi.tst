<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>LinSysIterDemo for a small linear problem of dimension 4</description>
    <executable>LinSysIterDemo</executable>
    <parameters>LinSysIterDemo_LinearSys_Jacobi.xml</parameters>
    <files>
        <file description="Session File">LinSysIterDemo_LinearSys_Jacobi.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">9.88991e-01</value>
            <value variable="v" tolerance="1e-12">2.01141e+00</value>
            <value variable="w" tolerance="1e-12">-1.01029e+00</value>
            <value variable="x" tolerance="1e-12">1.02135e+00</value>
        </metric>
    </metrics>
</test>
