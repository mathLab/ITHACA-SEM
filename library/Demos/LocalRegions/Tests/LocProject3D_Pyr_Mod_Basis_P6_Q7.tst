<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Project3D Pyramid Modified basis P=6 Q=7</description>
    <executable>LocProject3D</executable>
    <parameters>6 4 4 6 6 6 6 7 7 6  0 0 0  1 0 0  1 1 0  0 1 0  0.5 0.5 0.866</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-8">1.80317e-12</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-8">1.20881e-11</value>
        </metric>
    </metrics>
</test>
