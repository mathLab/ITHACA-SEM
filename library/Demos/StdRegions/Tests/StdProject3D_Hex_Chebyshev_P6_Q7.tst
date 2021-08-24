<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject3D Hexahedron Chebyshev basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s hexahedron -b chebyshev chebyshev chebyshev -o 6 6 6 -p 7 7 7</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">4.11256e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">4.1743e-13</value>
        </metric>
    </metrics>
</test>


