<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>LinSysIteratDemo for a small linear problem of dimension 4</description>
    <executable>LinSysIteratDemo</executable>
    <parameters>LinSysIteratDemo_LinearSys_Jacobi.xml</parameters>
    <files>
        <file description="Session File">LinSysIteratDemo_LinearSys_Jacobi.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">9.88991302e-01</value>
            <value variable="v" tolerance="1e-12">2.01141473e+00</value>
            <value variable="w" tolerance="1e-12">-1.01028590e+00</value>
            <value variable="x" tolerance="1e-12">1.02135051e+00</value>
        </metric>
    </metrics>
</test>
