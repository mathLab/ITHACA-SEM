<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdInterpBasis Tet Modified basis P=4 Q=5</description>
    <executable>StdInterpBasis</executable>
    <parameters> -s tetrahedron -b Modified_A Modified_B Modified_C -o 4 4 4 -p 5 5 5</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">1.89694e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">1.77636e-14</value>
        </metric>
    </metrics>
</test>
